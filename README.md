# TopoPerturb

Does knocking down a gene change the **shape** of a single-cell population — a loop
or a split into clusters — beyond a simple shift in mean or covariance?

TopoPerturb compares control and perturbed cells from CausalBench (K562, RPE1) with
**persistent homology on spectral graph distances** (Damrich, Berens & Kobak,
NeurIPS 2024), against conventional two-sample tests (energy distance). It is a
small, readable pipeline built around three habits:

1. **Get the data right first** — one verified source, an explicit audit, stated limits.
2. **Verify every method against a library**, and keep the three scientific questions
   (*differ? / beyond Gaussian? / what shape changed?*) separate.
3. **Calibrate on negative controls** before trusting any detection.

> A significant diagram difference is exploratory evidence of a structural change,
> not a validated biological loop. See [`docs/METHODS.md`](docs/METHODS.md).

**Current status.** Methods verified (8/8 checks). PCA is fit on a reserved half of
the controls that no test may draw. Calibration at 200 draws: all six methods hold
their 5% false-positive rate (mean 0.042). Synthetic: spectral H1 finds a ring that
energy misses. Real data: energy detects every K562 knockdown; RPE1 ZMAT2, CPSF3 and
the positive control have a robust `diffusion_H1` difference that a Gaussian
surrogate still explains. That is exploratory, not a validated loop.

## Layout

```text
topoperturb/        the package
  data.py           CausalBench load + PCA on a reserved reference-control subset
  distances.py      diffusion & effective-resistance distances on a kNN graph
  homology.py       persistence diagrams (ripser) + diagram distance (persim)
  twosample.py      energy / permutation / Gaussian-surrogate tests
  features.py       interpretable loop readout (lifetime, signed change, bootstrap CI)
  plotting.py       figures
experiments/        thin drivers; each saves one table + one figure to results/
  01_prepare_data   fetch, audit, 20/30/50-PC clouds
  02_synthetic      synthetic benchmark with known truth
  03_real           real K562 / RPE1 knockdowns
  04_calibration    control-vs-control false-positive rates
  05_compare_spaces 30-PC vs raw gene-space baseline
  verify_methods    independent reference checks
docs/               CONTEXT.md (living source of truth), DATA.md, METHODS.md,
                    UPDATE_REPORT.md (audit + calibration write-up)
results/            tables and figures (git-tracked)
data/               prepared clouds (gitignored; regenerate with step 1)
archive/            earlier exploratory phases (not imported by the pipeline)
```

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[verify]"
python experiments/verify_methods.py      # 8 checks vs networkx / definitions
```

Pinned versions live in `pyproject.toml` (Python 3.9).

## Reproduce

Each experiment is one command and writes its outputs to `results/`. The prepared
PCA clouds live in `data/` and are **gitignored**, so after a clone you must run
step 1 before anything else. The expression caches in `archive/data/` *are* in the
repo, so step 1 can stay offline.

Run them in this order (one at a time — each uses 8 worker processes):

```bash
# 1. Build the 30-PC clouds from the archived expression matrices (~seconds)
python experiments/01_prepare_data.py --cache 'archive/data/expression_{cell_line}.npz'

# 2c. Calibration — already done at 200 draws; rerun only if you change the clouds
python experiments/04_calibration.py          # ~20 min; -> results/calibration.*

# 2a. Synthetic benchmark with known truth (geometry vs topology)
python experiments/02_synthetic.py            # ~12 min; -> results/synthetic.{csv,png}

# 2b. Real K562 / RPE1 knockdowns + loop readout and diagrams
python experiments/03_real.py                 # ~8 min; -> results/real.{csv,png}, real_features.csv, diagrams/
```

Optional extras, same style:

```bash
# other PC counts (tagged files; do not overwrite the primary 30-PC clouds)
python experiments/01_prepare_data.py --n-pcs 20 --cache 'archive/data/expression_{cell_line}.npz'
python experiments/04_calibration.py --n-pcs 20

# 30-PC space vs raw gene space (same cells and permutations)
python experiments/05_compare_spaces.py --cache 'archive/data/expression_{cell_line}.npz'

# Hugging Face instead of the archive cache (downloads once, checksum-verified)
python experiments/01_prepare_data.py
```

Steps 2–4 take tens of minutes at full settings on 8 cores. Pass a smaller repeat
count (e.g. `python experiments/02_synthetic.py 10`) for a quick look, and
`--output results/local` to avoid overwriting saved results.

## Representation

PCA is **optional for the method** — every distance and test takes any point cloud —
and **required as a choice of input space**. The primary representation is 30 PCs
fitted on a reserved **reference half of the controls**; those cells are then
excluded from every comparison (`data.reference_split`). Fitting on the cells later
tested would make the null depend on the representation.

A 5-draw comparison against raw gene space (622 / 383 genes) keeps 30 PCs as
primary: about twice the distance contrast, about 1.4× faster, no null inflation.
Gene space is retained as a sensitivity baseline. Details: `docs/UPDATE_REPORT.md`.

## What to read first

- [`docs/CONTEXT.md`](docs/CONTEXT.md) — the living single-source-of-truth: data,
  statistics, method parameters, results, papers, and versions, with a changelog.
- [`docs/UPDATE_REPORT.md`](docs/UPDATE_REPORT.md) — PCA-independence audit,
  Damrich verification status, 20/30/50-PC comparison, calibration numbers.
- [`docs/DATA.md`](docs/DATA.md) — what the data is, the upstream efficacy filter,
  and why there is no guide/batch stratification.
- [`docs/METHODS.md`](docs/METHODS.md) — every formula, its reference, and how it is
  checked.
- `results/calibration.png` — control-vs-control false-positive rates (done).
- `results/synthetic.png` and `results/real.png` — headline figures, once steps 2a/2b
  have been run.

## Data

Hugging Face [`huggingbahl21/causalbench`](https://huggingface.co/datasets/huggingbahl21/causalbench),
pinned revision `068c7902f1740b8d653033c3decfed2c7156fa3d`. Prepared tables: already
library-normalised and `log1p`, efficacy-filtered. The only grouping label is
`intervention` — no guide, batch, or barcode, so permutations cannot be stratified.
Full limits in `docs/DATA.md`.
