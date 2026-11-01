# TopoPerturb — Running Context (single source of truth)

**Purpose.** One living document an agent (or the author) can read to understand the
whole project and make decisions: data, statistics, methods, results, papers, and
versions. Keep it current.

- **Last updated:** 2026-11-01
- **Repo state:** clean rewrite complete; experiments verified offline (8/8 method
  checks); PCA fit on a reserved reference half of the controls (§4); calibration run
  at 50 draws and PCA-vs-gene-space compared (§7a–b); **full-scale runs not yet
  executed**.
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

**Status: full-scale runs PENDING.** Only small verification/method-decision runs
executed so far. Full audit write-up: [`docs/UPDATE_REPORT.md`](UPDATE_REPORT.md).

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

### 7c. Earlier smoke observations

Observed in smoke runs (not final — low replicate counts, and these were run on the
**old balanced-fit PCA**, so even as smoke signals they must be re-run on the
controls-fit clouds before being read at all):
- Synthetic: `diffusion_H1`/`effres_H1` detect the ring where `energy` does not; null and
  pure translation stay low for PH.
- Real: `energy` rejects all knockdowns; some genes (e.g. RPL21, CPSF3, pos. control)
  also show `diffusion_H1`; control-vs-control stays ~0.
- Features: `signed_change` positive for the ring (loop appeared), as expected.

→ When full runs complete, paste the headline numbers here (detection rates, BH q-values,
calibration false-positive rates with Wilson CIs) and link the figures.

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
- [ ] Run full-scale `02/03` on the reference-fit clouds; fill §7 with final numbers
      + figures (the old smoke observations do not carry over).
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
