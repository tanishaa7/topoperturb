# TopoPerturb — Running Context (single source of truth)

**Purpose.** One living document an agent (or the author) can read to understand the
whole project and make decisions: data, statistics, methods, results, papers, and
versions. Keep it current.

- **Last updated:** 2026-11-01
- **Repo state:** methods verified (8/8); PCA on reserved reference controls (§4);
  **200-draw calibration, synthetic benchmark and 15-draw real knockdowns done**
  (§7a, §7c, §7d).
- **How to maintain:** when data, parameters, methods, results, or versions change,
  update the relevant section *and* add a dated line to the Changelog at the bottom.
  Numbers here must come from a real run — cite the file they came from.

---

## 1. Goal / question

Does knocking down a gene change the **shape** of a single-cell population (a loop /
H1 feature, or a split into clusters / H0) **beyond** a shift in mean or covariance?
Approach: persistent homology on spectral graph distances vs conventional two-sample
tests, on CausalBench K562/RPE1. Target: a clean, reproducible SciPy-style study repo.
A significant diagram difference is **exploratory** evidence, not a validated loop.

## 2. Repository layout

```
topoperturb/   package: data, distances, homology, twosample, features, plotting
experiments/   01_prepare_data, 02_synthetic, 03_real, 04_calibration,
               05_compare_spaces, verify_methods, config
docs/          DATA.md, METHODS.md, CONTEXT.md (this file), UPDATE_REPORT.md
results/       tables + figures (per experiment)
data/          cache/ (HF downloads), k562/, rpe1/ (prepared clouds) — all gitignored
archive/       earlier pipeline + recovered caches (NOT used by current code)
```

## 3. Data

- **Source:** Hugging Face `huggingbahl21/causalbench`, prepared tables (one TSV/line).
- **Pinned revision:** `068c7902f1740b8d653033c3decfed2c7156fa3d`
- **SHA-256:** k562 `455024ef…6fee5`, rpe1 `6ddd100e…cc408` (see `topoperturb/data.py`).
- **Preprocessing (upstream):** library-normalised + `log1p`; **efficacy-filtered**
  (strong perturbations; perturbed cells selected by target-gene response).
- **Only grouping label is `intervention`** → **no guide ID, no batch, no barcode.**
  Consequence: cannot stratify permutations by guide/batch; no per-guide validation
  from these inputs. (Full rationale in `docs/DATA.md`.)

### Selected cells (from `01_prepare_data`, 30-PC clouds)

| Cell line | cells | genes | controls (PCA reference / testing) | perturbations studied |
|---|---|---|---|---|
| K562 | 13,534 | 622 | 10,691 (5,346 / 5,345) | 8 pilot KDs |
| RPE1 | 13,899 | 383 | 11,485 (5,742 / 5,743) | 4 shared KDs + 1 positive control |

Half the controls are reserved to fit the PCA and never tested (§4); only the testing
half is available to steps 2b/2c. `results/data_summary.csv` carries both counts
(`n_cells`, `n_for_testing`).

**K562 per-group cells** (`results/data_summary.csv`): non-targeting 10,691; NCBP2 765;
CPSF3 575; ZMAT2 348; DLD 283; POLR2B 234; GFM1 230; POLRMT 208; RPL21 200.
**RPE1 per-group cells:** non-targeting 11,485; GFM1 1,669; POLRMT 212; ZMAT2 192;
CPSF3 173; pos. control (ENSG00000174243) 168.

**Gene IDs → symbols:** POLRMT `ENSG00000099821`, GFM1 `…168827`, ZMAT2 `…146007`,
CPSF3 `…119203`, RPL21 `…122026`, POLR2B `…047315`, DLD `…091140`, NCBP2 `…114503`,
pos. control `ENSG00000174243`.

## 4. Representation & control statistics

- **PCA:** 30 PCs, `fit="reference"`, seed 42, deterministic.
- **Reference split (PCA independence).** `data.reference_split` reserves **half the
  controls** to fit the axes and those cells are then excluded from every test;
  `data.testing_controls` defines the pool the tests may draw, and the mask travels
  with the clouds (`pca_reference`). K562: 5,346 reference / 5,345 testing controls.
  RPE1: 5,742 / 5,743. Draws are at most 300 cells, so the testing pool is ample.
- **Why not simply fit on all controls:** that removes perturbation leakage but
  leaves the null dependent on the representation — in a control-vs-control test the
  permuted cells are the cells that defined the coordinates, so its false-positive
  rate is not the one a fresh comparison would face. With a reserved half, nothing on
  either side of a comparison helped build the space. The split is identical across
  PC counts (it depends only on labels and the seed), so representations are compared
  on the same cells.
- **Earlier defaults, now sensitivity options only:** `"balanced"` was leakage-prone
  (~85–89% perturbation cells in the fit; its 30-D subspace overlapped the
  controls-only one by ~35%); `"controls"`/`"pooled"` sit between that and
  `"reference"`.
- Source: `data/<line>/clouds.json`, alternates `clouds_{20,50}pc.json` (step 1).

Primary representation (30 PCs; spectrum computed on the **testing** controls):

| Cell line | top-2 PC var frac | participation ratio | median centroid shift | var explained, 30 PCs |
|---|---|---|---|---|
| K562 | 0.280 | 15.00 | 1.78 | 15.1% |
| RPE1 | 0.330 | 12.49 | 3.00 | 20.9% |

**Representation sensitivity (20 / 30 / 50 PCs).** Prepared via
`01_prepare_data --n-pcs {20,50}` into tagged files:

| Cell line | PCs | var explained | top-2 var frac | participation ratio | median centroid shift |
|---|---|---|---|---|---|
| K562 | 20 | 12.0% | 0.338 | 10.78 | 1.74 |
| K562 | 30 | 15.1% | 0.280 | 15.00 | 1.78 |
| K562 | 50 | 20.9% | 0.211 | 24.42 | 1.86 |
| RPE1 | 20 | 16.3% | 0.416 | 8.31 | 2.98 |
| RPE1 | 30 | 20.9% | 0.330 | 12.49 | 3.00 |
| RPE1 | 50 | 28.9% | 0.238 | 21.94 | 3.05 |

The modest variance fraction is a property of sparse scRNA counts, not lost signal:
the control covariance spectrum flattens after ~PC 10 (K562 eigenvalues 2.60, 1.55,
0.80 at PC 5, 0.38 at PC 10, 0.28 at PC 20, 0.25 at PC 30, 0.20 at PC 50), so the
rest is thin, near-isotropic variance spread over hundreds of gene dimensions. Two
readings of the table support that: participation ratio grows roughly linearly with
the PC count (extra PCs add isotropic noise, not structure), while the median
perturbation centroid shift barely moves (1.74→1.86 K562, 2.98→3.05 RPE1), so the
perturbation signal is already captured by 20 PCs. 30 sits on a plateau.

(The K562 control covariance spectrum is what the synthetic benchmark uses as its
realistic base Gaussian.)

## 5. Methods & parameters

| Constant | Value | Meaning | Where |
|---|---|---|---|
| `N_PCS` | 30 | principal components | `data.embed` |
| `N_CELLS` | 150 | cells per group per comparison | `config` |
| `K_NN` | 15 | kNN graph neighbours (auto-raised if disconnected) | `distances` |
| `T_DIFF` | 8 | diffusion time | `distances` |
| `P_PERM` | 199 | permutation shuffles (p-resolution 1/200) | `twosample` |
| `M_SURR` | 99 | Gaussian surrogate pairs | `twosample` |
| `SEED` | 42 | global seed | everywhere |

**Method set (6):** `energy`, `covariance` (conventional); `euclidean_H0`,
`euclidean_H1`, `diffusion_H1`, `effres_H1` (persistent homology). Surrogate-adjusted
variants exist for each PH feature.

**Three separated questions** (`docs/METHODS.md`):
1. differ at all? → permutation test.
2. beyond a Gaussian? → covariance-matched surrogate test.
3. *what* shape changed? → `features.py`: leading H1 lifetime, **signed** change
   (appeared vs vanished), bootstrap 95% CI.

**Distances:** diffusion & effective-resistance on the kNN graph (Damrich et al. 2024);
effres uses the von Luxburg correction. Persistence via `ripser` (H0/H1), diagram
distance via `persim.wasserstein`. Distances scaled by mean off-diagonal so size ≠ shape.

**Numerical guards:** near-degenerate graphs/surrogates raise benign NumPy-2 FP flags;
guarded with `np.errstate` + `nan_to_num`, with a finiteness check on PCA output.

**Verification (`verify_methods.py`, independent references; all 8 checks pass):**
effective resistance vs `networkx.resistance_distance` → max error **1.1e-16**;
diffusion distance vs a direct `D⁻¹A` transition-matrix computation → **corr = 1.000000**
at t = 1/4/8 (the two forms differ by a global scale, which the mean-off-diagonal
scaling divides out); both spectral distances symmetric, zero-diagonal, non-negative;
a clean circle gives one dominant H1 loop (lifetime **1.32**); the PCA reference and
testing controls partition the control pool with no overlap and `embed`'s axes match a
PCA fitted directly on the reference rows; energy distance matches its definition;
diagram distance zero to itself, symmetric, monotone in bar length; permutation null
rejected **2/40 (~5%)**.

## 6. Experiments (what each produces)

| Script | Default scale | Outputs |
|---|---|---|
| `01_prepare_data` | both lines | `data/<line>/clouds.{npz,json}`, `results/data_summary.{csv,png}` |
| `02_synthetic` | 40 reps × 12 conditions | `results/synthetic.{csv,png}`, `synthetic_draws.csv` |
| `03_real` | 15 draws × 15 groups | `results/real.{csv,png}`, `real_draws.csv`, `real_features.csv`, `diagrams/` |
| `04_calibration` | 200 draws × {75,100,150} | `results/calibration.{csv,png}`, `calibration_draws.csv` |
| `05_compare_spaces` | 5 draws × 3 groups × 2 spaces | `results/space_comparison.csv`, `space_diagnostics.csv` |

Synthetic conditions: null, translation, scaling, anisotropy, heavy_tails, partial_kd
(mixture), ring, bimodal, ring_noise_{0.1,0.2,0.3,0.5}.

Note on `04_calibration`: its draws reuse cells from one control pool, so they are not
independent and the reported Wilson interval **understates** the true uncertainty — read
it as a lower bound. It draws only from the testing controls, so the PCA axes are
independent of the null being measured. `01_prepare_data --cache` accepts a
`{cell_line}` placeholder so both lines can be prepared in one run (the cross-line
summary needs that); `--n-pcs` on steps 1 and 2c selects an alternative representation
and tags its outputs (`clouds_50pc.npz`, `calibration_50pc.csv`, …).

## 7. Results

**Status:** calibration (200 draws), synthetic benchmark (40 reps) and real knockdowns
(15 draws) are **done**. Full audit write-up: [`docs/UPDATE_REPORT.md`](UPDATE_REPORT.md).

### 7a. Calibration on the revised representation (`results/calibration.csv`)

**Final run:** `04_calibration --draws 200` on the reference-fit 30-PC clouds, drawing
only from the testing controls; 1,200 draws, no missing p-values, 1,149 s. 36 rows
(2 lines × 3 sizes × {75,100,150} × 6 methods):

- false-positive rate at p < 0.05: **mean 0.042**, range 0.015–0.075;
- **0/36** rows with a Wilson interval lying entirely above 0.05 — no method is
  anti-conservative at any size in either line;
- 2/36 rows lie entirely below 0.05 (euclidean_H0 and euclidean_H1, RPE1, n=100, both
  0.015): slightly conservative, which costs power but not validity;
- per-method mean FPR: covariance 0.048, diffusion_H1 0.043, effres_H1 0.042,
  energy 0.044, euclidean_H0 0.042, euclidean_H1 0.032;
- worst row: covariance, RPE1, n=100, 0.075 (Wilson 0.046–0.120, includes 0.05);
- p-value uniformity: **1/36** KS tests below 0.05 (covariance, K562, n=75, p=0.034),
  under the ~1.8 expected by chance.

**Conclusion: all six methods are calibrated on the revised representation.** The
50-draw preliminary run's one outlier (energy, RPE1, n=75, 0.12) came down to 0.065
at 200 draws, i.e. it was sampling noise. Caveat that remains: the draws share one
control pool, so the Wilson intervals are lower bounds on the true uncertainty.

### 7b. Is PCA the right input space? (`results/space_{comparison,diagnostics}.csv`)

`05_compare_spaces --draws 5`: identical cells, permutations and methods in 30-PC
space vs raw gene space (622/383 genes). PCA is optional for the method — everything
after step 1 is dimension-agnostic — so this is a drop-in substitution.

| | 30 PCs | genes |
|---|---|---|
| distance concentration (std/mean of pairwise distances; higher = more contrast) | **0.291** / **0.224** | 0.136 / 0.108 |
| k needed to connect the kNN graph | 15 | 15 |
| seconds per draw (K562 / RPE1) | **12.8 / 10.8** | 20.1 / 14.3 |

(K562 / RPE1 in each cell.) The 15-NN sets agree only 24% (K562) and 40% (RPE1)
between the spaces, so PCA genuinely changes the graph the homology is built on.
Nulls stay ≤0.2 in both spaces; covariance is markedly more sensitive on K562 in PC
space (1.0 vs 0.2–0.4 of draws). One point against PCA: RPE1 POLRMT fires for energy
in 4/5 gene-space draws and 0/5 PC draws — at 5 draws, a hint, not a finding. The PH
features fire too rarely in both spaces at this scale to separate them.

**Decision: keep 30 PCs (reference-fit) as primary**, on distance contrast, runtime
and the absence of null inflation; gene space is retained as a sensitivity baseline.
Leading-H1 lifetimes differ a lot between spaces (K562 control 0.119 vs 0.011), so
always report the space alongside a loop readout.

### 7c. Synthetic benchmark (`results/synthetic.csv`, `results/synthetic.png`)

Full run: `02_synthetic.py` default, 40 replicates × 12 conditions, 480 draws, 1,085 s,
no missing p-values. Detection rate = fraction of replicates with permutation p < 0.05.

| Condition | energy | covariance | euclidean_H0 | euclidean_H1 | diffusion_H1 | effres_H1 |
|---|---|---|---|---|---|---|
| null | 0.05 | 0.02 | 0.00 | 0.00 | 0.08 | 0.05 |
| translation | **1.00** | 0.00 | 0.08 | 0.00 | 0.00 | 0.02 |
| scaling | **1.00** | **1.00** | 0.05 | 0.20 | 0.05 | 0.12 |
| anisotropy | 0.02 | 0.50 | **0.98** | 0.57 | 0.05 | 0.38 |
| heavy_tails | 0.45 | 0.05 | **1.00** | **1.00** | 0.10 | 0.85 |
| partial_kd | **1.00** | 0.12 | 0.05 | 0.00 | 0.05 | 0.00 |
| **ring** | 0.02 | 0.02 | 0.38 | 0.02 | **0.75** | **0.85** |
| bimodal | 0.02 | 0.02 | 0.10 | 0.05 | 0.30 | 0.18 |
| ring_noise_0.1 | 0.22 | 0.05 | 0.82 | 0.90 | **1.00** | **1.00** |
| ring_noise_0.2 | 0.08 | 0.05 | 0.30 | 0.02 | **1.00** | **1.00** |
| ring_noise_0.3 | 0.10 | 0.05 | 0.12 | 0.05 | **0.88** | 0.32 |
| ring_noise_0.5 | 0.08 | 0.05 | 0.02 | 0.05 | 0.08 | 0.05 |

**Headline.** Spectral H1 finds a ring that energy distance misses (0.75 / 0.85 vs
0.02), while ignoring a pure translation that energy always flags (0.00 vs 1.00).
That is the method's intended job. Diffusion H1 stays specific: it does not fire on
scaling, anisotropy, heavy tails or a partial knockdown, whereas Euclidean PH does.
Noise: both spectral methods are perfect at 0.1–0.2, diffusion still 0.88 at 0.3,
both collapse at 0.5. Null rates sit at the 5% they were calibrated to.

Caveat: `effres_H1` fires on heavy tails (0.85) and somewhat on anisotropy (0.38),
so it is less specific than diffusion. Prefer `diffusion_H1` as the primary shape
statistic; keep effective resistance as a corroborating spectral method.

### 7d. Real knockdowns (`results/real.csv`, `results/real.png`, `real_features.csv`)

Full run: `03_real.py` default, 15 draws × 15 groups, 225 tests, 666 s, 13 diagrams,
no missing p-values. Table = fraction of draws with permutation p < 0.05; primary
draw is BH-corrected across genes within a cell line (pos. control and
control-vs-control excluded from BH).

| Cell line | group | energy | covariance | euclidean_H0 | euclidean_H1 | diffusion_H1 | effres_H1 |
|---|---|---|---|---|---|---|---|
| K562 | POLRMT | **1.00** | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| K562 | GFM1 | **1.00** | 0.00 | 0.07 | 0.00 | 0.00 | 0.00 |
| K562 | ZMAT2 | **1.00** | 0.33 | 0.00 | 0.00 | 0.00 | 0.00 |
| K562 | CPSF3 | **1.00** | 0.33 | 0.00 | 0.00 | 0.13 | 0.13 |
| K562 | RPL21 | **1.00** | **0.93** | 0.73 | 0.47 | 0.13 | 0.07 |
| K562 | POLR2B | **1.00** | **1.00** | 0.73 | 0.20 | 0.27 | 0.67 |
| K562 | DLD | **1.00** | 0.47 | 0.00 | 0.00 | 0.13 | 0.00 |
| K562 | NCBP2 | **1.00** | **0.87** | 0.33 | 0.20 | 0.40 | 0.00 |
| K562 | control_vs_control | 0.00 | 0.13 | 0.13 | 0.07 | 0.07 | 0.07 |
| RPE1 | ZMAT2 | **1.00** | **0.87** | 0.73 | 0.60 | **0.93** | 0.67 |
| RPE1 | CPSF3 | **1.00** | **1.00** | 0.67 | 0.33 | **1.00** | **0.80** |
| RPE1 | POLRMT | 0.33 | 0.13 | 0.00 | 0.13 | 0.00 | 0.07 |
| RPE1 | GFM1 | 0.33 | 0.07 | 0.13 | 0.13 | 0.00 | 0.07 |
| RPE1 | pos. control | **1.00** | **1.00** | 0.40 | 0.60 | **1.00** | 0.67 |
| RPE1 | control_vs_control | 0.00 | 0.00 | 0.00 | 0.07 | 0.00 | 0.00 |

BH q on the primary draw, genes only, q < 0.05:

- **energy:** every K562 gene (q = 0.005); RPE1 ZMAT2, CPSF3 (q = 0.010) and GFM1
  (q = 0.020). RPE1 POLRMT is not (p = 0.070).
- **covariance:** K562 RPL21, POLR2B (q = 0.020); RPE1 ZMAT2, CPSF3 (q = 0.020).
- **diffusion_H1:** K562 RPL21 (q = 0.040, but only 13% of draws — not robust);
  RPE1 ZMAT2 and CPSF3 (q = 0.010, 93% / 100% of draws). Pos. control also 100% of
  draws (excluded from BH).
- **effres_H1:** RPE1 ZMAT2 (q = 0.020). POLR2B is 67% of draws but primary q = 0.120.

**Reading, three questions kept separate.**

1. *Differ at all?* Yes, for every K562 knockdown and for RPE1 ZMAT2, CPSF3, GFM1 and
   the positive control. Energy is the workhorse; control-vs-control is 0.00.
2. *Beyond a Gaussian?* The covariance-matched surrogate **never** rejects
   `diffusion_H1` (0/15 draws on every group, including the RPE1 hits). A Gaussian
   with the same mean and covariance can produce the same diagram distances, so a
   diagram difference is **not** evidence of a non-Gaussian loop.
3. *What shape?* Leading-H1 signed changes are small and their bootstrap intervals
   all cover both signs. Do not claim a loop appeared or vanished.

So the RPE1 ZMAT2 / CPSF3 / pos. control `diffusion_H1` hits are **repeatable
differences in the persistence diagram**, calibrated and robust across draws, but
they sit inside a Gaussian-matched null. They are exploratory structure-beyond-mean
signals, not validated loops. K562 spectral hits are mostly not robust across draws;
POLR2B / RPL21 / NCBP2 look like covariance / stretch (energy + covariance fire,
spectral H1 does not stably).

## 8. Papers / references

- **Damrich, Berens & Kobak (2024)** "Persistent homology for high-dimensional data
  based on spectral methods," *NeurIPS 2024* — diffusion & effective-resistance PH (core).
- **Székely & Rizzo** — energy distance (conventional comparator).
- **von Luxburg et al.** — effective-resistance correction on large graphs.
- **CausalBench (2025)**, *Commun. Biol.* — dataset + upstream efficacy filter.
- **scPerturb (2024)**, *Nat. Methods* — energy statistics for perturbation comparison.
- **SCEPTRE / Barry et al. (2024)**, *Genome Biol.* — calibration in single-cell CRISPR.
- **Fasy et al. (2014)** — confidence sets for persistence diagrams.
- Future comparators (not implemented): Pertpy, scGeom, CONCORD.

## 9. Software versions

Python **3.9.6**; numpy **2.0.2**, scipy **1.13.1**, scikit-learn **1.6.1**,
pandas **2.3.3**, matplotlib **3.9.4**, ripser **0.6.15**, persim **0.4.1**,
huggingface_hub (latest); networkx (verify-only). Pinned in `pyproject.toml`.

## 10. Open decisions / TODO

- [x] `04_calibration --draws 200` at 30 PCs — done, all methods calibrated (§7a).
- [ ] `04_calibration --draws 200 --n-pcs 20` and `--n-pcs 50` for stability.
- [x] `02_synthetic` at 40 reps — done (§7c). Spectral H1 finds the ring; energy does not.
- [x] `03_real 15` — done (§7d). Energy hits all K562 KDs; RPE1 ZMAT2/CPSF3/pos. control
      have robust `diffusion_H1`, but the Gaussian surrogate does not reject.
- [ ] `05_compare_spaces --draws 25` — is the RPE1 POLRMT gene-space/PCA gap real?
- [ ] Close the Damrich et al. gap: side-by-side against the authors' reference
      implementation, and re-check their `t`/`k` defaults and input space against the
      paper (our checks are against independent definitions, not their code).
- [ ] Report the `fit="balanced"` (and `drop_genes`) runs as an explicit representation
      sensitivity check alongside the reference-fit primary analysis.
- [ ] Fetch canonical data from Hugging Face (step 1 without `--cache`) and confirm
      clouds match the offline-cache clouds.
- [ ] Decide whether to add a modest novel method (currently: referenced methods only).
- [ ] Interpret top real hits biologically (cell-cycle/stress programs) before any
      loop claim; current evidence is exploratory.
- [ ] If stratified/robust calibration is needed, obtain guide/batch metadata upstream.

## 11. Changelog

- **2026-11-01** — Clean rewrite: package + experiments; legacy moved to `archive/`;
  methods verified against references; docs (DATA/METHODS/CONTEXT) added. Full-scale
  runs pending.
- **2026-11-01** — Audit + fixes. (1) **PCA default `balanced` → `controls`**: the
  balanced fit was ~85–89% perturbation cells (controls 11% K562 / 17% RPE1) and its
  30-D subspace overlapped the controls-only one by only ~35%, i.e. the axes were
  partly defined by the cells under test (double dipping); clouds regenerated for both
  lines, §4 stats and `results/data_summary.{csv,png}` updated. (2) `verify_methods.py`
  extended from 4 to 7 checks (diffusion vs transition matrix, distance-matrix
  properties, circle → one dominant H1; diagram check strengthened to symmetry +
  monotonicity) — all pass. (3) `twosample.surrogate_test` guards the z-score against a
  zero-variance surrogate spread. (4) Calibration's non-independent draws documented as
  making its Wilson interval a lower bound. (5) `01_prepare_data --cache` takes a
  `{cell_line}` placeholder so both lines (and thus the cross-line summary) come from
  one run. Full-scale runs still pending; §7 smoke observations predate the PCA change.
- **2026-11-01** — PCA independence, representation audit, first calibration
  (`docs/UPDATE_REPORT.md`). (1) **`fit="controls"` → `fit="reference"`**: half the
  controls are reserved to fit the axes and excluded from all testing
  (`data.reference_split`, `data.testing_controls`, `pca_reference` mask in the
  clouds), so no cell both defines and is judged in the space; verified by an 8th
  check in `verify_methods.py`. §4 stats regenerated. (2) **20/30/50-PC clouds**
  prepared (`--n-pcs`, tagged outputs); 30 PCs sits on a plateau — spectrum flattens
  after ~PC 10, participation ratio grows linearly, centroid shift unchanged. (3)
  **Calibration run** (50 draws) on the revised representation: mean FPR 0.041,
  uniformity failures at chance level — §7a. (4) **New `05_compare_spaces`**: 30 PCs
  vs raw gene space on identical cells/permutations; PCA kept as primary on distance
  contrast (~2×), runtime (~1.4×) and no null inflation — §7b. Still pending: 200-draw
  calibration, full `02`/`03`, and a side-by-side against the Damrich et al. code.
- **2026-11-01** — **200-draw calibration complete** (30 PCs, 1,149 s): mean FPR 0.042,
  0/36 rows significantly above 5%, 1/36 KS uniformity failures. All six methods
  calibrated; the 50-draw outlier was noise. §7a and `UPDATE_REPORT.md` §6 updated.
- **2026-11-01** — README rewritten: it still contained an unresolved merge conflict
  between the archived poster pipeline (`topoperturb_benchmark/`, `earlier_phases/`)
  and the current package. Stale `requirements.txt` (old pins + kmapper) now points
  at `pyproject.toml`; `CLAUDE.md` lists `05_compare_spaces`.
- **2026-11-01** — **Synthetic benchmark complete** (40 reps, 1,085 s): spectral H1
  finds the ring (diffusion 0.75, effres 0.85) that energy misses (0.02) and ignores
  a pure translation (0.00 vs energy 1.00). Diffusion is more specific than
  effective resistance (no fire on heavy tails). §7c. `03_real` still pending.
- **2026-11-01** — **Real knockdowns complete** (15 draws, 666 s): energy q < 0.05 on
  every K562 gene and RPE1 ZMAT2/CPSF3/GFM1; `diffusion_H1` robust on RPE1 ZMAT2,
  CPSF3 and pos. control (93–100% of draws) but the Gaussian surrogate never
  rejects. Control-vs-control energy is 0.00. §7d.
