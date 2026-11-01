# Data: what it is, and what it is not

The analysis uses the **prepared** CausalBench tables published on Hugging Face
(`huggingbahl21/causalbench`, one gzipped TSV per cell line), pinned to a single
revision and verified by SHA-256 in [`topoperturb/data.py`](../topoperturb/data.py).

## What one table contains

- One **cell per row**; one **gene per column**, plus a final `intervention` column.
- `intervention` is the knocked-down gene (an Ensembl ID) or `non-targeting` for a
  control. A small number of rows are labelled `excluded` and are dropped.
- Values are already **library-normalised and `log1p`-transformed** counts
  (finite, non-negative). We add no normalisation, log, imputation or filtering.

`data.audit()` scans the whole table once and refuses it unless every value is
finite and non-negative, every cell has a label, and controls are present.

## Upstream selection (important for interpretation)

The published tables are **already efficacy-filtered**: CausalBench selected strong
perturbations and kept perturbed cells partly by their target-gene response. A
detected difference therefore describes *selected responders*, not every cell
assigned a guide. The control-vs-control calibration (step 4) does **not**
reproduce this selection mechanism, so it bounds only the ordinary false-positive
rate, not any bias introduced by responder selection.

## What the table does **not** contain

- **No guide ID, no batch, no cell barcode.** The only grouping is `intervention`.
- `source_rows` saved with the clouds is the zero-based row index in the TSV, kept
  for traceability — it is *not* a biological identifier.

Consequence: we cannot stratify permutations by guide or batch, and cannot do
per-guide validation, from these inputs alone. Doing so would need richer upstream
data. This is a limitation we state, not one we have worked around.

## From cells to clouds (step 1)

`experiments/01_prepare_data.py` keeps the controls plus the perturbations we study
(the eight K562 pilot knockdowns; for RPE1 the four shared knockdowns plus one
strong positive control) and projects every cell onto **30 principal components**.

PCA is fit on a **held-out reference half of the controls** (`fit="reference"`).
`data.reference_split` splits the control pool in two with a fixed seed: one half
defines the axes and is then **retired from the analysis**, the other half is the
only control material any test may draw. `data.testing_controls` is the single
place that pool is defined, and `01_prepare_data` stores the mask in the clouds so
steps 2b and 2c exclude exactly the cells the axes were fitted on.

Why go further than fitting on all controls. Fitting on every control removes
*perturbation* leakage but not the subtler kind: in a control-vs-control null the
permuted cells are the same cells that defined the coordinates, so the null
distribution is not independent of the representation being tested, and the
false-positive rate it reports is not quite the one a fresh comparison would see.
With a reserved reference half, nothing on either side of a comparison helped build
the space. The cost is cells, and we can afford it: ~5.3k testing controls in K562
and ~5.7k in RPE1 against draws of at most 300.

Earlier defaults, kept only as sensitivity options: `fit="balanced"` was outright
leakage-prone (with many small perturbation groups, equal-per-group sampling made
the fit ~85–89% perturbation cells, and its 30-D subspace overlapped the
controls-only one by only ~35%); `fit="controls"` and `fit="pooled"` sit between
that and the reference split. We fix one primary representation deliberately rather
than searching over PCA/graph settings for favourable p-values; `data.embed` exposes
`n_pcs`, `fit` and `drop_genes` for anyone who wants to probe sensitivity, and
`01_prepare_data --n-pcs` writes alternative representations to tagged files so a
sensitivity run cannot overwrite the primary clouds. The control covariance
spectrum saved alongside the clouds — computed on the **testing** controls, since
those are the cells the tests actually see — is what the synthetic benchmark (step
2a) uses to make its base Gaussian realistic.

### Why 30 PCs, when they explain only ~15–21% of variance

The low fraction is expected for sparse single-cell counts and is not a sign that
structure was left behind: the control covariance spectrum flattens after roughly
the tenth component. In K562 the eigenvalues run 2.60, 1.55, 0.80 (PC 5), 0.38
(PC 10), 0.28 (PC 20), 0.25 (PC 30), 0.20 (PC 50) — PCs 20–50 are a near-flat
noise floor, so the residual variance is spread thinly over hundreds of gene
dimensions rather than concentrated in a few we are ignoring. Two consequences
visible in the prepared clouds: the participation ratio grows almost linearly with
the PC count (K562 10.8 → 15.0 → 24.4 at 20/30/50 PCs), i.e. extra components add
near-isotropic noise rather than new structure; and the median perturbation
centroid shift is essentially unchanged (K562 1.74 / 1.78 / 1.86; RPE1 2.98 / 3.00
/ 3.05), i.e. the perturbation signal sits in the leading components and is already
captured at 20. Thirty is therefore a point on a plateau rather than a tuned
choice, and `--n-pcs 20` / `--n-pcs 50` exist to re-check that claim.
