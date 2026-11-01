"""Load the CausalBench perturbation data and turn it into PCA point clouds.

The inputs are the *prepared* CausalBench tables published on Hugging Face: one
cell line per gzipped TSV, already library-normalised and ``log1p`` transformed,
with the upstream efficacy filter applied (strong perturbations, and perturbed
cells selected by their target-gene response). Each row is one cell; every column
is a gene except the last, ``intervention``, which names the knocked-down gene or
``non-targeting`` for a control. No normalisation, imputation or extra filtering
is done here -- see ``docs/DATA.md`` for what that means for interpretation.

Typical use (see ``experiments/01_prepare_data.py``)::

    path = fetch("k562")                 # download + checksum, or reuse the cache
    genes, counts, report = audit(path)
    X, labels, rows = select_cells(path, targets, min_cells=150)
    Z, pca = embed(X, labels, n_pcs=30)  # one label-blind PCA per cell line
"""
import gzip
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

# The one pinned source. A checksum guards against a silently changed table.
REPO = "huggingbahl21/causalbench"
REVISION = "068c7902f1740b8d653033c3decfed2c7156fa3d"
SHA256 = {
    "k562": "455024efc9f40b877f7f9c1d55f9eaaf957c7995c6289ffbcd65251fc856fee5",
    "rpe1": "6ddd100eff9dbcab65ff94b662a87e9eb50de6401b5242c0f0108692c92cc408",
}
CONTROL = "non-targeting"
CACHE = Path(__file__).resolve().parents[1] / "data" / "cache"


def sha256(path):
    """Return the hex SHA-256 of a file, read in large blocks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(cell_line, cache=CACHE):
    """Return a local path to the prepared TSV, downloading it once if needed.

    The download is verified against the pinned checksum before it is accepted,
    so a corrupt or wrong-revision file never reaches the analysis.
    """
    cache.mkdir(parents=True, exist_ok=True)
    name = f"causalbench_{cell_line}.tsv.gz"
    path = cache / name
    if not path.exists():
        url = f"https://huggingface.co/datasets/{REPO}/resolve/{REVISION}/data/prepared/{name}"
        partial = path.with_suffix(path.suffix + ".part")
        with urlopen(url, timeout=120) as response, partial.open("wb") as out:
            shutil.copyfileobj(response, out)
        if sha256(partial) != SHA256[cell_line]:
            partial.unlink()
            raise ValueError(f"Downloaded {name} does not match the pinned checksum")
        partial.replace(path)
    if sha256(path) != SHA256[cell_line]:
        raise ValueError(f"Cached {name} does not match the pinned checksum")
    return path


def gene_columns(path):
    """Read the header and return the gene names (every column but the last)."""
    with gzip.open(path, "rt") as handle:
        names = handle.readline().rstrip("\r\n").split("\t")
    if len(names) < 2 or names[-1] != "intervention" or len(set(names)) != len(names):
        raise ValueError("Expected unique gene columns followed by 'intervention'")
    return names[:-1]


def _chunks(path, genes):
    """Stream the TSV in row blocks so a multi-hundred-MB table fits in memory."""
    types = {gene: np.float32 for gene in genes}
    types["intervention"] = "string"
    return pd.read_csv(path, sep="\t", compression="gzip", dtype=types, chunksize=10_000)


def audit(path):
    """Scan the whole table once and confirm it is what we expect.

    Checks every value is finite and non-negative (prepared log1p counts), every
    cell carries a label, and controls are present. Returns the gene list, a
    per-intervention cell count, and a small summary so the data step is explicit.
    """
    genes, counts = gene_columns(path), Counter()
    low, high, zero_rows = np.inf, -np.inf, 0
    for chunk in _chunks(path, genes):
        if chunk.intervention.isna().any() or chunk.intervention.str.strip().eq("").any():
            raise ValueError("Missing intervention label")
        X = chunk[genes].to_numpy()
        if not np.isfinite(X).all() or (X < 0).any():
            raise ValueError("Prepared expression must be finite and non-negative")
        counts.update(chunk.intervention.tolist())
        low, high = min(low, float(X.min())), max(high, float(X.max()))
        zero_rows += int((X.sum(axis=1) == 0).sum())
    if not counts or counts[CONTROL] == 0:
        raise ValueError("Empty table or no non-targeting controls")
    report = {"n_cells": sum(counts.values()), "n_genes": len(genes),
              "expression_min": low, "expression_max": high, "all_zero_cells": zero_rows}
    return genes, counts, report


def eligible_targets(counts, min_cells):
    """Perturbations with enough cells to compare against controls."""
    return sorted(g for g, n in counts.items()
                  if g not in {CONTROL, "excluded"} and n >= min_cells)


def select_cells(path, targets, min_cells=150, min_controls=300):
    """Load controls plus the requested perturbations into one array.

    ``targets`` are Ensembl IDs. Returns the expression matrix, the matching
    intervention labels, and the zero-based source rows (so a cell can always be
    traced back to the table -- note the table has no cell barcode or batch).
    """
    genes, counts, _ = audit(path)
    if counts[CONTROL] < min_controls:
        raise ValueError(f"Need {min_controls} controls; found {counts[CONTROL]}")
    targets = list(dict.fromkeys(targets))
    for gene in targets:
        if gene in {CONTROL, "excluded"} or counts.get(gene, 0) < min_cells:
            raise ValueError(f"Target {gene}: needs >= {min_cells} cells; found {counts.get(gene, 0)}")
    keep = set(targets) | {CONTROL}
    blocks, labels, rows = [], [], []
    for chunk in _chunks(path, genes):
        mask = chunk.intervention.isin(keep)
        blocks.append(chunk.loc[mask, genes].to_numpy(np.float32))
        labels.extend(chunk.loc[mask, "intervention"].tolist())
        rows.extend(chunk.index[mask].tolist())
    X = np.concatenate(blocks)
    return X, np.asarray(labels), np.asarray(rows, dtype=np.int64)


def load_cache_npz(path, keep=None):
    """Fast path: load an already-parsed ``(X, lab, genes)`` expression cache.

    Lets the pipeline run offline from a saved matrix instead of re-reading the
    large TSV; the PCA and all downstream steps are identical either way.
    """
    with np.load(path, allow_pickle=False) as z:
        X, labels, genes = z["X"], z["lab"].astype(str), z["genes"].astype(str)
    mask = labels != "excluded"
    if keep is not None:
        mask &= np.isin(labels, list(keep))
    return X[mask], labels[mask], list(genes)


def reference_split(labels, frac=0.5, seed=42, min_testing=300):
    """Reserve a disjoint subset of controls to fit the PCA on, and nothing else.

    Fitting the axes on controls that later enter a test leaves the test's own null
    distribution dependent on the representation being tested: the permuted labels
    are shuffled within cells that helped define the coordinates. Holding out a
    reference subset removes that dependence -- the axes come from control cells no
    test ever sees, so a control-vs-control comparison is a genuine null.

    Returns a boolean mask over ``labels`` (``True`` = reference, PCA-fit only).
    """
    control_idx = np.flatnonzero(labels == CONTROL)
    n_reference = int(round(frac * len(control_idx)))
    if len(control_idx) - n_reference < min_testing:
        raise ValueError(f"Reserving {n_reference} of {len(control_idx)} controls for the PCA "
                         f"leaves fewer than {min_testing} for testing")
    rng = np.random.default_rng([seed, 101])
    mask = np.zeros(len(labels), dtype=bool)
    mask[rng.choice(control_idx, n_reference, replace=False)] = True
    return mask


def testing_controls(labels, reference):
    """Indices of the controls the tests may draw from: those not reserved for the PCA."""
    return np.flatnonzero((labels == CONTROL) & ~np.asarray(reference, dtype=bool))


def embed(X, labels, n_pcs=30, fit="reference", reference=None, drop_genes=None,
          genes=None, seed=42):
    """Project cells onto ``n_pcs`` principal components.

    ``fit`` chooses which cells define the axes. The default ``"reference"`` fits on
    the held-out reference controls alone (see :func:`reference_split`), so neither
    the perturbation under test nor the control cells it is compared against had any
    hand in building the representation. ``"controls"`` fits on all controls -- no
    perturbation leakage, but the test controls do shape their own axes, so it is a
    sensitivity option rather than the primary choice. ``"pooled"`` fits on every
    cell (controls dominate by count, so it is close to ``"controls"``);
    ``"balanced"`` fits on an equal number of cells per group, which -- with many
    small perturbation groups -- makes the fit mostly perturbation cells and is kept
    only for sensitivity checks. ``drop_genes`` optionally removes target genes
    before PCA to test whether a result survives without the direct knockdown signal.
    """
    rng = np.random.default_rng(seed)
    if drop_genes:
        feature = ~np.isin(genes, list(drop_genes))
        X = X[:, feature]
    if fit == "pooled":
        fit_idx = np.arange(len(X))
    elif fit == "reference":
        if reference is None:
            reference = reference_split(labels, seed=seed)
        fit_idx = np.flatnonzero(reference)
    elif fit == "controls":
        fit_idx = np.flatnonzero(labels == CONTROL)
    elif fit == "balanced":
        groups = [np.flatnonzero(labels == g) for g in sorted(set(labels))]
        n = min(len(g) for g in groups)
        fit_idx = np.sort(np.concatenate([rng.choice(g, n, replace=False) for g in groups]))
    else:
        raise ValueError("fit must be 'reference', 'controls', 'balanced' or 'pooled'")
    variable = X[fit_idx].var(axis=0) > 0        # drop constant genes so PCA is well posed
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed)
    # Randomized SVD on float32 raises benign floating-point flags under NumPy 2;
    # the finiteness check below is the real guard on the output.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        pca.fit(X[np.ix_(fit_idx, np.flatnonzero(variable))])
        Z = pca.transform(X[:, variable]).astype(np.float64)
    if not np.isfinite(Z).all():
        raise ValueError("PCA produced non-finite coordinates")
    return Z, pca


def control_spectrum(Z, labels, reference=None):
    """Describe the control cloud: the summary a synthetic benchmark needs.

    Returns the control covariance eigenvalues, its participation ratio (an
    anisotropy measure), the fraction of variance in the top two PCs, and the
    median centroid shift of the perturbations -- a realistic translation size.
    Pass ``reference`` to describe the *testing* controls only, which is what the
    tests actually see (the PCA-fit controls are reserved and never compared).
    """
    control = Z[labels == CONTROL] if reference is None else Z[testing_controls(labels, reference)]
    cov = np.cov(control.T)
    eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]
    shifts = [np.linalg.norm(Z[labels == g].mean(0) - control.mean(0))
              for g in sorted(set(labels)) if g != CONTROL]
    pr = float(eigenvalues.sum() ** 2 / (eigenvalues ** 2).sum())
    return {"eigenvalues": eigenvalues.tolist(), "participation_ratio": pr,
            "top2_var_frac": float(eigenvalues[:2].sum() / eigenvalues.sum()),
            "median_centroid_shift": float(np.median(shifts)) if shifts else 0.0}


def save_clouds(path, Z, labels, rows, meta, reference=None):
    """Write the PCA clouds and their provenance/summary side by side.

    ``reference`` is the PCA-fit control mask, saved with the clouds so every
    downstream step excludes exactly the cells the axes were fitted on.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if reference is None:
        reference = np.zeros(len(Z), dtype=bool)
    np.savez_compressed(path, Z=Z, labels=labels, source_rows=rows,
                        pca_reference=np.asarray(reference, dtype=bool))
    path.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")


def load_clouds(path):
    """Load ``(Z, labels, source_rows, pca_reference), meta`` from :func:`save_clouds`."""
    path = Path(path)
    with np.load(path, allow_pickle=False) as z:
        clouds = z["Z"], z["labels"].astype(str), z["source_rows"], z["pca_reference"]
    meta = json.loads(path.with_suffix(".json").read_text())
    return clouds, meta
