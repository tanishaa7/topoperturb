# CLAUDE.md

**Read [`docs/CONTEXT.md`](docs/CONTEXT.md) first.** It is the single source of truth
for this project: data (source, revision, counts), statistics, method parameters,
results, papers, and software versions. Keep it current — when data, parameters,
methods, results, or versions change, update the relevant section and add a dated
Changelog line.

Quick orientation:
- Package: `topoperturb/` (data, distances, homology, twosample, features, plotting).
- Experiments: `experiments/01_prepare_data … 04_calibration`, plus `verify_methods`.
- Method formulas/references and the three separated questions: `docs/METHODS.md`.
- Data provenance and its limits (no guide/batch): `docs/DATA.md`.
- `archive/` is earlier work, not imported by the current pipeline.

Environment: `.venv/` (Python 3.9, pinned stack in `pyproject.toml`).
