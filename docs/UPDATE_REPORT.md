# Update report — PCA independence, representation choice, calibration

**Date:** 2026-11-01. Scope: audit the existing spectral-PH pipeline, fix the
representation's independence from the tests, check whether PCA is the right input
space at all, and get calibration on the record before any findings are read.
No new algorithms were added. All numbers below name the file they came from.

**Status of the numbers:** the calibration (§6) is the full 200-draw run and can be
quoted. The space comparison (§4) is 5 draws per group — a small, reproducible run
meant to support a *method* decision, not a finding. The real-data and synthetic
results are still pending.

---

## 1. Audit: where PCA is required, and where it is not

Tracing expression → clouds → tests:

| Stage | Code | Needs PCA? |
|---|---|---|
| select cells, log1p (upstream), efficacy filter | `data.select_cells` | no |
| project to 30 PCs | `data.embed` | **this *is* the PCA step** |
| kNN graph | `distances.knn_graph` | no — takes any point cloud |
| diffusion / effective-resistance distance | `distances.*` | no — operates on the graph |
| persistence diagrams | `homology.diagrams` | no — takes a distance matrix |
| permutation / surrogate tests | `twosample.*` | no — takes two clouds |
| synthetic base Gaussian | `data.control_spectrum` | yes, by construction (it describes the PCA cloud) |

So PCA is **optional for the method and load-bearing only as a choice of input
space**. Everything after step 1 is dimension-agnostic, which is what makes the
gene-space baseline in §3 a drop-in substitution rather than a parallel pipeline.
The one real coupling is the synthetic benchmark (step 2a), which models the control
covariance of whatever representation it is handed.

## 2. Fix: PCA independence via a reserved reference-control subset

**Problem.** Fitting the axes on controls that later enter a test leaves the test's
null dependent on the representation being tested: in a control-vs-control
comparison the permuted cells are the same cells that defined the coordinates, so
the measured false-positive rate is not the one a fresh comparison would face. (The
earlier `fit="balanced"` default was worse still — ~85–89% of the fit cells were
perturbation cells, and its 30-D subspace overlapped the controls-only one by only
~35%.)

**Fix.** `data.reference_split` reserves half the controls to fit the PCA and
retires them from the analysis; `data.testing_controls` is the only definition of
the pool the tests may draw; the mask ships with the clouds (`pca_reference`) so
steps 2b/2c exclude exactly the fit cells. K562 5,346 reference / 5,345 testing,
RPE1 5,742 / 5,743, against draws of at most 300 cells — the cost in cells is
affordable. The split depends only on labels and the seed, so it is identical across
PC counts and representations stay comparable on the same cells.

**Verified** in `verify_methods.py` (`check_pca_reference_split`): the reference and
testing sets partition the control pool with empty intersection, no perturbed cell
ever enters the fit, the split is deterministic, and `embed`'s axes equal a PCA
fitted directly on the reference rows alone.

## 3. Verification status against Damrich et al. (2024)

| Component | Reference checked against | Result |
|---|---|---|
| effective resistance (Laplacian pseudoinverse) | `networkx.resistance_distance` | max error 1.1e-16 |
| diffusion distance | direct `D⁻¹A` transition-matrix computation, t = 1, 4, 8 | corr = 1.000000 (the two forms differ by a global scale, divided out by the mean-off-diagonal scaling) |
| distance-matrix validity | symmetry, zero diagonal, non-negativity | pass |
| H1 detection | clean circle → one dominant loop | lifetime 1.32, next bar >5× smaller |
| diagram distance | self-distance, symmetry, monotonicity in bar length | pass |

**What is implemented:** the two spectral distances (diffusion with the eigen-form,
effective resistance with the von Luxburg degree correction), Vietoris–Rips
persistence via `ripser`, Wasserstein diagram distance via `persim`, mean-off-diagonal
scaling so size differences do not read as shape differences.

**Limitation, stated plainly:** these are checks against *independent definitions*,
not against the authors' own code. We have not run the paper's reference
implementation side by side, and the paper's own parameter choices (notably whether
it applies these distances to PCA-reduced or full-dimensional single-cell input, and
its `t` and `k` defaults) have not been re-checked against the text in this pass.
That comparison is the honest remaining gap in "verified against Damrich et al."

## 4. Comparison: 30 PCs vs raw gene space

`experiments/05_compare_spaces.py` runs the identical battery, the identical cells and
the identical permutations in both spaces (`results/space_comparison.csv`,
`results/space_diagnostics.csv`). Gene space is the normalised log1p matrix as
published — 622 genes (K562), 383 (RPE1).

Graph and geometry, one 150-cell control draw:

| Cell line | space | dims | k to connect | mean degree | distance concentration | 15-NN agreement with the other space |
|---|---|---|---|---|---|---|
| K562 | 30 PCs | 30 | 15 | 22.8 | **0.291** | 0.24 |
| K562 | genes | 622 | 15 | 27.8 | 0.136 | 0.24 |
| RPE1 | 30 PCs | 30 | 15 | 22.7 | **0.224** | 0.40 |
| RPE1 | genes | 383 | 15 | 26.9 | 0.108 | 0.40 |

Distance concentration is `std/mean` of the pairwise distances: **higher is better**
(more contrast between near and far). Gene space roughly halves it, which is the
high-dimensional concentration effect that motivates spectral distances in the first
place. Both spaces connect at k = 15, so neither is pathological. The 15-NN agreement
of only 24% (K562) and 40% (RPE1) shows PCA is **not** a cosmetic step: it changes
which cells count as neighbours, and therefore the graph the homology is built on.

Battery over 5 paired draws (fraction with p < 0.05; control-vs-control rows are
false positives, knockdown rows sensitivity):

| Cell line | group | method | genes | 30 PCs |
|---|---|---|---|---|
| K562 | CPSF3 | covariance | 0.2 | **1.0** |
| K562 | NCBP2 | covariance | 0.4 | **1.0** |
| K562 | CPSF3 / NCBP2 | energy | 1.0 | 1.0 |
| K562 | control_vs_control | all | ≤0.2 | ≤0.2 |
| RPE1 | POLRMT | energy | 0.8 | 0.0 |
| RPE1 | GFM1 | energy | 1.0 | 0.8 |
| RPE1 | control_vs_control | all | ≤0.2 | 0.0 |

Runtime per draw (mean): K562 12.8 s in 30 PCs vs 20.1 s in gene space; RPE1 10.8 s
vs 14.3 s. Leading diffusion-H1 lifetimes differ substantially between spaces
(e.g. K562 control 0.119 in PCs vs 0.011 in genes), i.e. the loop readout is
representation-dependent and must always be reported with its space.

**Reading.** PCA earns its place on distance contrast (≈2× better), runtime (≈1.4×
faster) and the absence of any null inflation, and it does not cost sensitivity for
the conventional comparators on K562. The one point *against* it is RPE1 POLRMT,
where energy fires in 4/5 gene-space draws and 0/5 PCA draws — a hint that some
perturbation signal lives outside the leading 30 PCs. At 5 draws that is a single
observation, not a result. **No claim is made that either space is superior for the
PH features**: at this draw count the PH methods fire too rarely in both spaces to
separate them.

## 5. Representation robustness: 20 / 30 / 50 PCs

Prepared with `01_prepare_data --n-pcs {20,50}` into tagged files; see `CONTEXT.md`
§4 for the table. The short version: cumulative explained variance is 12.0 / 15.1 /
20.9% (K562) and 16.3 / 20.9 / 28.9% (RPE1), which looks low but is a property of
sparse scRNA counts — the control covariance spectrum flattens after ~PC 10 (K562
eigenvalues 2.60, 1.55, 0.80 at PC 5, 0.38 at PC 10, 0.28 at PC 20, 0.25 at PC 30,
0.20 at PC 50). Two checks say 30 sits on a plateau rather than at a tuned optimum:
the participation ratio grows almost linearly with the PC count (10.8 → 15.0 → 24.4
for K562), i.e. added components contribute near-isotropic noise; and the median
perturbation centroid shift barely moves (1.74 → 1.86 K562, 2.98 → 3.05 RPE1), i.e.
the perturbation signal is already in the leading components.

**Not yet done:** calibration and the battery at 20 and 50 PCs. The clouds exist, so
this is one command per PC count (§7).

## 6. Calibration on the revised representation

`04_calibration --draws 200` on the reference-fit 30-PC clouds, drawing only from the
testing controls (`results/calibration.csv`; 1,200 draws, 36 cell-line × size × method
rows, no missing p-values):

- False-positive rate at p < 0.05: **mean 0.042**, range 0.015–0.075.
- **0 of 36** rows have a Wilson interval lying entirely above 0.05.
- 2 of 36 lie entirely below it (euclidean_H0/H1, RPE1, n=100, both 0.015) —
  slightly conservative, a power cost rather than a validity problem.
- Per method, mean FPR: covariance 0.048, diffusion_H1 0.043, effres_H1 0.042,
  energy 0.044, euclidean_H0 0.042, euclidean_H1 0.032.
- p-value uniformity (KS vs Uniform(0,1)): **1 of 36** below 0.05, under the ~1.8
  chance predicts.

**Reading.** All six methods hold their nominal 5% rate on the revised
representation, including both spectral PH features. A preliminary 50-draw run had
one outlier (energy, RPE1, n=75, 0.12); at 200 draws it is 0.065, so it was sampling
noise. One caveat stands: the draws reuse cells from one control pool, so the Wilson
intervals **understate** uncertainty and should be read as lower bounds.

## 7. Recommendation and what remains

**Recommended primary method:** the existing pipeline, unchanged in substance —
30 PCs fitted on reserved reference controls, kNN graph at k = 15, diffusion and
effective-resistance persistent homology, permutation tests for significance with the
covariance-matched surrogate as the "beyond Gaussian" check, energy distance as the
conventional comparator. Report gene space as a sensitivity baseline, not as the
primary, and always state the space alongside any H1 lifetime.

Remaining work, in order:

1. ~~`04_calibration --draws 200`~~ — done (§6).
2. `04_calibration --draws 200 --n-pcs 20` and `--n-pcs 50` — calibration stability
   across representations.
3. `03_real` and `02_synthetic` at full scale on the reference-fit clouds; §7 of
   `CONTEXT.md` is still empty and the old smoke observations predate this change.
4. `05_compare_spaces --draws 25` — enough draws to say whether the RPE1 POLRMT
   gene-space/PCA discrepancy in §4 is real.
5. Side-by-side check against the Damrich et al. reference implementation, closing
   the gap named in §3.
