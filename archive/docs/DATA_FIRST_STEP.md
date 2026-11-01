# First step: inspect, select, then choose a representation

This is the new **prepared-table route**. `src/prepare_data.py` remains the legacy
benchmark reconstruction. It still uses the old caches and should not be used to
evaluate these new options. No archived data/results are replaced.

## What we verified about the source

Repository: https://huggingface.co/datasets/huggingbahl21/causalbench

Pinned revision: `068c7902f1740b8d653033c3decfed2c7156fa3d`.
The downloader checks the SHA-256 advertised for each large file at this revision.
See its [dataset card](https://huggingface.co/datasets/huggingbahl21/causalbench/blob/068c7902f1740b8d653033c3decfed2c7156fa3d/README.md)
and [upstream code](https://huggingface.co/datasets/huggingbahl21/causalbench/blob/068c7902f1740b8d653033c3decfed2c7156fa3d/upstream/preprocessing.py).

The prepared TSVs already incorporate:

1. Selection of strong perturbations: >50 DEGs, knockdown <= -0.3,
   and >25 filtered cells in the summary workbook.
2. Selection of cells whose measured target expression is <= the control
   10th percentile (for targets that can be measured).
3. Removal of zero-total cells, per-cell normalization to the retained-cell
   median total using all original genes, and natural `log1p`.
4. Restriction to measured genes targeted in >100 retained cells. Unsupported
   target labels become `excluded`; those cells remain in the exported table.

These are source-code/card observations, not a new raw-data reproduction here.
The tables have no original cell barcodes, guide IDs or batch metadata. The
source TSV row is a reproducible identifier, but cannot recover lost metadata.
The downloaded file inventory contains TSVs, not the NPZ outputs mentioned as
rebuild products in the card. Raw H5AD files total about 19 GB.

**Do not normalize or log these values again.** Their row totals cover a restricted
gene panel, not all original transcripts. Re-normalizing them would change the
scientific quantity. Nor can this route undo upstream efficacy selection.

## Run the small scripts

Use the project's Python environment with `requirements.txt` installed. The
downloader uses Python's standard library; no Hugging Face SDK is required.

```bash
# 1. Download one verified table (RPE1 ~57 MB; K562 ~119 MB).
python src/fetch_inputs.py rpe1

# 2. Audit all rows without fitting PCA or picking targets.
python src/prepare_inputs.py data/cache/causalbench/causalbench_rpe1.tsv.gz \
  --cell-line rpe1 --audit-only --output results/local/rpe1_audit

# 3. Read groups.csv. Explicitly choose the four historical shared targets:
# POLRMT, GFM1, ZMAT2, CPSF3. This is a continuity example, not a new discovery set.
python src/prepare_inputs.py data/cache/causalbench/causalbench_rpe1.tsv.gz \
  --cell-line rpe1 --targets ENSG00000099821 ENSG00000168827 ENSG00000146007 ENSG00000119203 \
  --pca-fit balanced --n-pcs 30 --output results/local/rpe1_balanced
```

Use `k562` and its downloaded filename for K562. Choose IDs from that audit's
`groups.csv`; the scripts never assume the eight historical targets. For a broad
count-selected analysis, replace `--targets ...` with `--all-eligible`.
Defaults require >=150 cells per target and >=300 controls, supporting later
150-vs-150 comparisons and a disjoint control-vs-control comparison. This is
availability selection **within the already efficacy-selected source**.

Each command needs a new output directory; no overwrite is allowed. The loader
rejects files not matching the pinned source and chosen cell line. A different
dataset revision should be reviewed and its hashes updated deliberately.

## Choices to compare and why

| Choice | Option | Reason / limitation |
|---|---|---|
| Default expression | No extra transform | Already normalized and log-transformed upstream. |
| Balanced PCA | `--pca-fit balanced` | Equal fitting-cell counts per group prevent abundant controls dominating the basis. Default cap 150/group; all selected cells are then transformed. Still depends on the chosen groups. |
| Pooled PCA | `--pca-fit pooled` | Uses all selected cells; closest to the old benchmark. Group abundance affects the basis. |
| Control PCA | `--pca-fit controls` | Basis is fitted only on controls; may miss a perturbation-specific direction. |
| Gene scaling | `--gene-scaling none` (default) or `unit-variance` | Scaling equalizes feature variance but can amplify weak/noisy genes. Scale is estimated only on PCA fitting cells. PCA centers in both cases; it does not whiten PCs. |
| Target removal | `--drop-target-genes` | Drops every selected target gene from the shared feature set before PCA, testing dependence on direct knockdown. It cannot undo upstream cell selection. |
| PCA dimension | `--n-pcs 10`, `30`, `50` | Tests sensitivity to retained variance/noise. No claim that 30 is optimal. |
| Group selection | Explicit IDs or `--all-eligible` | No centroid-shift screen or automatically selected positive control. |

Start with balanced/pooled/control PCA at 30 PCs, no scaling, identical targets.
Then compare target removal against the baseline. Change one choice at a time,
and keep the seed/selected cells fixed. Avoid choosing a setting because it makes
the desired gene significant. This stage creates **no independent test split**;
it supports exploration and comparison, not confirmatory inference.

PCA drops genes constant in its fitting cells and reports them. The audit reports
zero rows in the restricted gene panel without silently discarding them: zero
in this panel is not evidence of a zero-total raw cell.

## What flows out

```text
verified prepared TSV: cells × genes + intervention
  -> audit all cells/labels/values
  -> keep explicit eligible targets + non-targeting controls
  -> optional target-gene removal, fitted feature scaling
  -> fit PCA on balanced / pooled / control cells
  -> transform every selected cell in the same basis
```

| Output | Contents |
|---|---|
| `groups.csv` | Counts for every source label, eligibility and selection flags; includes excluded rows in the audit. |
| `manifest.json` | Source revision/hash, options, selection, PCA fitting counts, explained variance and limitations. |
| `expression.npz` | Unmodified selected expression `X`, labels `lab`, original panel `genes`, source row indices. Target removal affects PCA, not this audit copy. |
| `clouds.npz` | `Z` (cells × PCs), `lab`, `source_rows`. One cell line per run. |
| `pca.npz` | Feature names, components, mean, scaling, explained variance and exact fitting row IDs. |

Audit-only writes the first two files. The one-line-per-run cloud format is
deliberately separate from the old `Zk/Zr` combined benchmark file. `real.py`
still reads `data/clouds.npz`; running it will **not** use this new output.
To inspect the input that a later comparison will receive:

```python
import numpy as np
z = np.load("results/local/rpe1_balanced/clouds.npz", allow_pickle=False)
A = z["Z"][z["lab"] == "non-targeting"]
B = z["Z"][z["lab"] == "ENSG00000146007"]  # ZMAT2
print(A.shape, B.shape)  # all available cells; methods.battery needs sampled groups
```

Before connecting to confirmatory tests, decide the discovery/evaluation split,
experimental blocking and multiple-testing family. Cell subsampling is not
biological replication. A new raw-count route is needed to compare efficacy
filter on/off, retain original metadata or choose a broader gene panel. That
route has not been implemented or run in this prepared-first step.
