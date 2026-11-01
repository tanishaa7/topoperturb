"""Audit prepared CausalBench TSVs, then select cells and fit an explicit PCA.

These inputs are already normalized/log1p transformed and efficacy-filtered.
No normalization, log transform, imputation or additional response filtering is
performed here. Output is one cell line per directory, separate from legacy data.
"""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from fetch_inputs import REPO, REVISION, SHA256, sha256

CONTROL = "non-targeting"


def columns(path):
    with gzip.open(path, "rt") as handle:
        names = handle.readline().rstrip("\r\n").split("\t")
    if (len(names) < 2 or names[-1] != "intervention"
            or len(set(names)) != len(names) or any(not name for name in names)):
        raise ValueError("Expected unique gene columns followed by intervention")
    return names[:-1]


def chunks(path, genes):
    types = {gene: np.float32 for gene in genes}
    types["intervention"] = "string"
    return pd.read_csv(path, sep="\t", compression="gzip", dtype=types, chunksize=10000)


def audit(path):
    genes, counts = columns(path), Counter()
    minimum, maximum, zero_rows = float("inf"), float("-inf"), 0
    for chunk in chunks(path, genes):
        if chunk.intervention.isna().any() or chunk.intervention.str.strip().eq("").any():
            raise ValueError("Missing intervention label")
        X = chunk[genes].to_numpy()
        if not np.isfinite(X).all() or (X < 0).any():
            raise ValueError("Prepared expression must be finite and nonnegative")
        counts.update(chunk.intervention.tolist())
        minimum, maximum = min(minimum, float(X.min())), max(maximum, float(X.max()))
        zero_rows += int((X.sum(axis=1) == 0).sum())
    if not counts or counts[CONTROL] == 0:
        raise ValueError("Empty table or no non-targeting controls")
    unknown = set(counts) - set(genes) - {CONTROL, "excluded"}
    if unknown:
        raise ValueError(f"Unsupported labels (not measured gene IDs): {sorted(unknown)}")
    return genes, counts, {"n_rows": sum(counts.values()), "n_genes": len(genes),
                           "minimum": minimum, "maximum": maximum,
                           "zero_expression_rows": zero_rows}


def select_targets(counts, targets, all_eligible, min_cells, min_controls):
    if counts[CONTROL] < min_controls:
        raise ValueError(f"Need {min_controls} controls; found {counts[CONTROL]}")
    if all_eligible:
        targets = sorted(g for g, n in counts.items()
                         if g not in {CONTROL, "excluded"} and n >= min_cells)
    targets = list(dict.fromkeys(targets or []))
    if not targets:
        raise ValueError("Choose --targets GENE_ID ... or --all-eligible; inspect groups.csv first")
    for gene in targets:
        if gene in {CONTROL, "excluded"} or counts[gene] < min_cells:
            raise ValueError(f"Target {gene}: requires >= {min_cells} cells; found {counts[gene]}")
    return targets


def load_selected(path, genes, targets):
    matrices, labels, source_rows = [], [], []
    keep = set(targets) | {CONTROL}
    for chunk in chunks(path, genes):
        mask = chunk.intervention.isin(keep)
        matrices.append(chunk.loc[mask, genes].to_numpy(np.float32))
        labels.extend(chunk.loc[mask, "intervention"].tolist())
        source_rows.extend(chunk.index[mask].tolist())
    return np.concatenate(matrices), np.asarray(labels), np.asarray(source_rows, dtype=np.int64)


def fit_pca(X, labels, mode, n_pcs, fit_cells, scaling, seed):
    rng = np.random.default_rng(seed)
    if mode == "pooled":
        fit = np.arange(len(X))
    elif mode == "controls":
        fit = np.flatnonzero(labels == CONTROL)
    else:
        groups = [np.flatnonzero(labels == g) for g in sorted(set(labels))]
        n = min(fit_cells, min(map(len, groups)))
        fit = np.sort(np.concatenate([rng.choice(g, n, replace=False) for g in groups]))
    variance = X[fit].var(axis=0, dtype=np.float64)
    retained = variance > 0
    if not retained.any():
        raise ValueError("No variable genes in the PCA fitting cells")
    if not 1 <= n_pcs <= min(len(fit) - 1, int(retained.sum())):
        raise ValueError("n-pcs exceeds the fitting-cell or variable-gene limit")
    scale = np.sqrt(variance[retained]) if scaling == "unit-variance" else np.ones(retained.sum())
    Y = X[:, retained] / scale
    model = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed)
    model.fit(Y[fit])
    Z = model.transform(Y)
    if not np.isfinite(Z).all():
        raise ValueError("PCA produced nonfinite coordinates")
    return Z, model, scale, fit, retained


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Prepared .tsv.gz (not raw counts)")
    parser.add_argument("--cell-line", choices=["k562", "rpe1"], required=True)
    parser.add_argument("--output", type=Path, required=True, help="New directory; refuses overwrite")
    parser.add_argument("--audit-only", action="store_true")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--targets", nargs="+", help="Explicit Ensembl IDs, not gene symbols")
    selection.add_argument("--all-eligible", action="store_true")
    parser.add_argument("--min-cells", type=int, default=150)
    parser.add_argument("--min-controls", type=int, default=300)
    parser.add_argument("--pca-fit", choices=["balanced", "pooled", "controls"], default="balanced")
    parser.add_argument("--fit-cells-per-group", type=int, default=150)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--gene-scaling", choices=["none", "unit-variance"], default="none")
    parser.add_argument("--drop-target-genes", action="store_true",
                        help="Remove ALL selected target genes from every cell before PCA")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if min(args.min_cells, args.min_controls, args.fit_cells_per_group, args.n_pcs) < 1 or args.seed < 0:
        parser.error("Counts/components must be positive; seed must be nonnegative")
    if args.output.exists():
        parser.error("Output already exists; choose a new directory")
    digest = sha256(args.input)
    # Do not silently attribute an arbitrary or wrong-cell-line table to this source.
    if digest != SHA256[args.cell_line]:
        parser.error("Input checksum differs from the pinned prepared table for this cell line")
    genes, counts, report = audit(args.input)
    targets = [] if args.audit_only else select_targets(
        counts, args.targets, args.all_eligible, args.min_cells, args.min_controls)
    inventory = pd.DataFrame([
        {"intervention": g, "n_cells": n,
         "eligible": g not in {CONTROL, "excluded"} and n >= args.min_cells,
         "selected": g in targets or (bool(targets) and g == CONTROL)}
        for g, n in sorted(counts.items())])
    manifest = {"repository": REPO, "revision": REVISION,
                "input_path": str(args.input.resolve()), "input_sha256": digest,
                "cell_line": args.cell_line, "audit": report,
                "options": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "selected_targets": targets, "normalization": "none added; prepared log1p values",
                "upstream_selection": "strong perturbations + target-expression efficacy filter",
                "row_identity": "zero-based source TSV row, NOT original cell barcode",
                "inference_scope": "exploratory; no independent evaluation split created"}
    if not args.audit_only:
        X, labels, rows = load_selected(args.input, genes, targets)
        feature_mask = ~np.isin(genes, targets) if args.drop_target_genes else np.ones(len(genes), bool)
        feature_genes = np.asarray(genes)[feature_mask]
        Z, pca, scale, fit, variable = fit_pca(
            X[:, feature_mask], labels, args.pca_fit, args.n_pcs,
            args.fit_cells_per_group, args.gene_scaling, args.seed)
        manifest.update({"n_selected_cells": len(X), "n_pca_features": int(variable.sum()),
                         "pca_fit_group_counts": dict(Counter(labels[fit].tolist())),
                         "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
                         "constant_genes_removed": feature_genes[~variable].tolist()})
    args.output.mkdir(parents=True)
    inventory.to_csv(args.output / "groups.csv", index=False)
    if not args.audit_only:
        np.savez_compressed(args.output / "expression.npz", X=X, lab=labels,
                            genes=np.asarray(genes), source_rows=rows)
        np.savez_compressed(args.output / "clouds.npz", Z=Z, lab=labels, source_rows=rows)
        np.savez_compressed(args.output / "pca.npz", genes=feature_genes[variable],
                            components=pca.components_, mean=pca.mean_, scale=scale,
                            explained_variance=pca.explained_variance_, fit_source_rows=rows[fit])
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), **report, "selected_targets": targets}, indent=2))


if __name__ == "__main__":
    main()
