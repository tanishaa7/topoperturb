# Archive — earlier work, not the current pipeline

These files are kept for provenance and are **not** part of the TopoPerturb
pipeline described in the top-level `README.md`. Nothing here is imported by the
`topoperturb` package or the `experiments/` scripts.

- `src/` — the earlier flat script pipeline and its `earlier/` phase-by-phase
  exploratory analyses. The current package reimplements the same methods cleanly.
- `tests/` — replay and formula checks tied to the earlier pipeline and recovered
  outputs.
- `data/` — recovered, pre-parsed expression matrices and PCA clouds, plus the
  original metadata/provenance. The current pipeline instead fetches the prepared
  tables from Hugging Face (see `docs/DATA.md`). The expression matrices here can
  still be used offline via `01_prepare_data.py --cache archive/data/expression_<line>.npz`.
- `docs/` — the earlier research status, data guide, verification record, and the
  historical `archive.zip`.
