"""Step 2d -- is the PCA step earning its place? Same tests, two input spaces.

Nothing downstream of step 1 needs PCA: ``distances`` and ``homology`` only ever see
a point cloud and return a distance matrix, so the gene-expression matrix can be fed
in directly. This script runs the *identical* battery, cells and permutations in both
spaces -- 30 PCs (the primary representation) and the raw normalised gene space (622
genes for K562, 383 for RPE1) -- so the only thing that differs is the input space.

Three things are reported:

* graph/geometry diagnostics: the k needed to connect the kNN graph, distance
  concentration (std/mean of pairwise distances -- lower means less contrast between
  near and far, the usual high-dimensional failure mode), and how much the 15-nearest
  neighbours of a cell agree between the two spaces;
* the test battery on a control-vs-control null (false positives) and on the two
  largest knockdowns per line (sensitivity), with the leading H1 lifetimes;
* wall-clock seconds per draw in each space.

This is a small, paired comparison, not a full-scale run: read it as a decision aid
about the representation, not as the study's findings.

    python experiments/05_compare_spaces.py --cache 'archive/data/expression_{cell_line}.npz'
    python experiments/05_compare_spaces.py --draws 10
"""
import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
from sklearn.neighbors import NearestNeighbors, kneighbors_graph

import config
from topoperturb import data, distances, features, twosample

CONTROL = data.CONTROL
N, SEED = config.N_CELLS, config.SEED
CELL_LINES = ("k562", "rpe1")
SPACES = ("pca30", "genes")
_CELLS = {}


def cells(cell_line, cache):
    """``(X_genes, Z_pca, labels, testing_control_idx)`` for one line, cached per process.

    The gene matrix and the PCA clouds are loaded through the same selection used in
    step 1, so row *i* is the same cell in both; the assertion below keeps that honest.
    """
    if cell_line not in _CELLS:
        targets = config.TARGETS[cell_line]
        keep = set(targets) | {CONTROL}
        path = cache.format(cell_line=cell_line) if cache else None
        if path:
            X, labels, _genes = data.load_cache_npz(path, keep=keep)
        else:
            X, labels, _rows = data.select_cells(data.fetch(cell_line), targets)
        (Z, cloud_labels, _rows, reference), _meta = data.load_clouds(config.cloud_path(cell_line))
        assert len(X) == len(Z) and (labels == cloud_labels).all(), "cell order must match"
        _CELLS[cell_line] = (X.astype(np.float64), Z, labels,
                             data.testing_controls(labels, reference))
    return _CELLS[cell_line]


def cloud(space, X, Z, idx):
    """The point cloud for one space, restricted to the drawn cells."""
    return Z[idx] if space == "pca30" else X[idx]


def draw_indices(cell_line, group, draw, control_idx, labels):
    """Cell indices for one (group, draw); identical in both spaces, so the test is paired.

    The seed is built from the group's position in ``config.TARGETS`` rather than from
    ``hash(group)``, which varies between worker processes and would silently unpair
    the two spaces.
    """
    order = ["control_vs_control"] + list(config.TARGETS[cell_line])
    rng = np.random.default_rng([SEED, 205, CELL_LINES.index(cell_line), order.index(group), draw])
    if group == "control_vs_control":
        pick = rng.choice(control_idx, 2 * N, replace=False)
        return pick[:N], pick[N:], rng
    return (rng.choice(control_idx, N, replace=False),
            rng.choice(np.flatnonzero(labels == group), N, replace=False), rng)


def run(task):
    """One (line, group, draw, space): the battery, the leading H1 change, and runtime."""
    cell_line, group, draw, space, cache = task
    X, Z, labels, control_idx = cells(cell_line, cache)
    a_idx, b_idx, rng = draw_indices(cell_line, group, draw, control_idx, labels)
    A, B = cloud(space, X, Z, a_idx), cloud(space, X, Z, b_idx)
    start = time.time()
    row = {"cell_line": cell_line, "group": config.symbol(group), "draw": draw,
           "space": space, "n_dims": A.shape[1]}
    for name, (p, obs) in twosample.permutation_test(
            twosample.all_statistics, A, B, rng=rng).items():
        row[f"p_{name}"] = p
    row.update(features.signed_change(A, B))
    row["seconds"] = time.time() - start
    return row


def _k_to_connect(P, k=distances.K_NN):
    """Smallest k (stepping by 5, as ``distances.knn_graph`` does) giving one component.

    This is graph quality in the sense that matters here: the spectral distances are
    undefined across components, so a space that needs a larger k is one where the
    neighbourhoods are less trustworthy.
    """
    while True:
        A = kneighbors_graph(P, k, mode="connectivity", include_self=False).toarray()
        A = np.maximum(A, A.T)
        if connected_components(A, directed=False)[0] == 1:
            return k, A
        k += 5


def diagnostics(cell_line, cache):
    """Graph and geometry diagnostics for one control draw, in both spaces."""
    X, Z, labels, control_idx = cells(cell_line, cache)
    a_idx, _b_idx, _rng = draw_indices(cell_line, "control_vs_control", 0, control_idx, labels)
    rows, neighbours = [], {}
    for space in SPACES:
        P = cloud(space, X, Z, a_idx)
        k_connected, A = _k_to_connect(P)
        off = pdist(P)
        neighbours[space] = NearestNeighbors(n_neighbors=distances.K_NN + 1).fit(P) \
            .kneighbors(P, return_distance=False)[:, 1:]
        rows.append({"cell_line": cell_line, "space": space, "n_dims": P.shape[1],
                     "k_to_connect": k_connected, "mean_degree": float(A.sum(1).mean()),
                     "distance_concentration": float(off.std() / off.mean())})
    shared = np.mean([len(set(a) & set(b)) / distances.K_NN
                      for a, b in zip(neighbours["pca30"], neighbours["genes"])])
    for row in rows:
        row["knn_agreement_pca_vs_genes"] = float(shared)
    return rows


def main(draws, cache, out):
    out.mkdir(parents=True, exist_ok=True)
    groups = {}
    for cell_line in CELL_LINES:
        _X, _Z, labels, _ctl = cells(cell_line, cache)
        counts = pd.Series(labels).value_counts()
        biggest = [g for g in counts.index if g != CONTROL][:2]      # the two best-powered KDs
        groups[cell_line] = ["control_vs_control"] + biggest

    diag = pd.DataFrame([r for cl in CELL_LINES for r in diagnostics(cl, cache)])
    diag.to_csv(out / "space_diagnostics.csv", index=False)
    print("\ngraph / geometry diagnostics (one 150-cell control draw):")
    print(diag.round(3).to_string(index=False))

    tasks = [(cl, g, d, s, cache) for cl in CELL_LINES for g in groups[cl]
             for d in range(draws) for s in SPACES]
    with Pool(8) as pool:
        table = pd.DataFrame(pool.map(run, tasks, chunksize=1))
    table.to_csv(out / "space_comparison.csv", index=False)

    detect = table.pivot_table(index=["cell_line", "group"], columns="space",
                               values=[f"p_{m}" for m in twosample.METHODS],
                               aggfunc=lambda p: (p < 0.05).mean())
    print(f"\nfraction of {draws} draws with p < 0.05 (control_vs_control rows are "
          "false positives, the rest sensitivity):")
    print(detect.round(2).to_string())
    print("\nleading diffusion-H1 lifetimes and runtime per draw:")
    print(table.groupby(["cell_line", "space"])[
        ["leading_H1_control", "leading_H1_perturbed", "leading_H1_change", "seconds"]]
        .mean().round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--draws", type=int, default=5, help="draws per group per space")
    parser.add_argument("--cache", type=str, help="parsed expression .npz ({cell_line} allowed)")
    parser.add_argument("--output", type=Path, default=config.RESULTS)
    args = parser.parse_args()
    t = time.time()
    main(args.draws, args.cache, args.output)
    print(f"{time.time() - t:.0f}s")
