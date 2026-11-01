"""Step 1 -- get the data right before anything else depends on it.

Fetches a CausalBench cell line from Hugging Face (checksum-verified), keeps the
controls and the perturbations we study, and projects every cell onto 30 PCs with
one PCA fitted on a held-out reference half of the controls -- cells that no later
test is allowed to draw, so the representation is independent of every comparison
made in it. Saves the clouds (with the reference mask) and a provenance/summary
record per line, plus a cross-line summary table and figure.

    python experiments/01_prepare_data.py                 # both lines, from Hugging Face
    python experiments/01_prepare_data.py --cell-line k562
    python experiments/01_prepare_data.py --cache archive/data/expression_k562.npz --cell-line k562
    python experiments/01_prepare_data.py --cache 'archive/data/expression_{cell_line}.npz'

The ``--cache`` option reads an already-parsed expression matrix instead of the
large TSV, so the pipeline runs offline and gives identical PCA clouds. A
``{cell_line}`` placeholder in the path is filled in per line, so both lines can
be prepared in one run (needed for the cross-line summary table and figure).
"""
import argparse
import json

import numpy as np
import pandas as pd

import config
from topoperturb import data, plotting


def get_cells(cell_line, targets, cache):
    """Return (expression, labels, source_rows, genes) for controls + targets."""
    keep = set(targets) | {data.CONTROL}
    if cache:
        X, labels, genes = data.load_cache_npz(cache, keep=keep)
        rows = np.arange(len(X))                       # the cache does not carry source rows
    else:
        path = data.fetch(cell_line)
        X, labels, rows = data.select_cells(path, targets)
        genes = data.gene_columns(path)
    return X, labels, rows, np.asarray(genes)


def prepare(cell_line, cache, n_pcs):
    """Build and save the PCA clouds for one cell line; return a summary row list."""
    cache = cache.format(cell_line=cell_line) if cache else None
    targets = config.TARGETS[cell_line]
    X, labels, rows, genes = get_cells(cell_line, targets, cache)
    reference = data.reference_split(labels, seed=config.SEED)
    Z, pca = data.embed(X, labels, n_pcs=n_pcs, fit="reference", reference=reference,
                        genes=genes, seed=config.SEED)
    spectrum = data.control_spectrum(Z, labels, reference=reference)
    meta = {"cell_line": cell_line, "source": config.data_source(cache),
            "n_cells": len(Z), "n_genes": int(genes.size),
            "targets": {g: config.symbol(g) for g in targets},
            "pca_fit": "reference", "n_pcs": n_pcs,
            "n_reference_controls": int(reference.sum()),
            "n_testing_controls": int(data.testing_controls(labels, reference).size),
            "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
            "control_spectrum": spectrum,
            "provenance": "prepared CausalBench (log1p, efficacy-filtered); labels only, no batch"}
    data.save_clouds(config.cloud_path(cell_line, n_pcs), Z, labels, rows, meta,
                     reference=reference)
    counts = pd.Series(labels).value_counts()
    held_out = {data.CONTROL: int(reference.sum())}
    return [{"cell_line": cell_line, "group": config.symbol(g), "n_cells": int(counts[g]),
             "n_for_testing": int(counts[g]) - held_out.get(g, 0)}
            for g in [data.CONTROL] + targets]


def main(cell_lines, cache, n_pcs):
    rows, specs = [], {}
    for cell_line in cell_lines:
        rows += prepare(cell_line, cache, n_pcs)
        (_Z, _lab, _rows, _ref), meta = data.load_clouds(config.cloud_path(cell_line, n_pcs))
        specs[cell_line] = meta
        print(f"{cell_line}: {meta['n_cells']} cells, {meta['n_pcs']} PCs, "
              f"{meta['n_reference_controls']} controls reserved for the PCA fit, "
              f"{meta['n_testing_controls']} left for testing; "
              f"top-2 PC variance {meta['control_spectrum']['top2_var_frac']:.2f}")
    summary = pd.DataFrame(rows)
    config.RESULTS.mkdir(parents=True, exist_ok=True)
    tag = "" if n_pcs == config.N_PCS else f"_{n_pcs}pc"
    summary.to_csv(config.RESULTS / f"data_summary{tag}.csv", index=False)

    fig, axes = plotting.plt.subplots(1, len(specs), figsize=(4.5 * len(specs), 3.2), squeeze=False)
    for ax, (cell_line, meta) in zip(axes[0], specs.items()):
        ratio = np.array(meta["explained_variance_ratio"])
        ax.bar(range(1, len(ratio) + 1), ratio * 100, color="#4477aa")
        ax.set_title(f"{cell_line}: PCA variance explained "
                     f"({100 * ratio.sum():.1f}% in {meta['n_pcs']} PCs)", fontsize=9)
        ax.set_xlabel("principal component")
        ax.set_ylabel("variance (%)")
    plotting.save(fig, config.RESULTS / f"data_summary{tag}.png")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cell-line", choices=["k562", "rpe1"], help="default: both")
    parser.add_argument("--cache", type=str, help="parsed expression .npz for offline use")
    parser.add_argument("--n-pcs", type=int, default=config.N_PCS,
                        help=f"principal components (default {config.N_PCS}; "
                             "other counts are written to tagged files)")
    args = parser.parse_args()
    main([args.cell_line] if args.cell_line else ["k562", "rpe1"], args.cache, args.n_pcs)
