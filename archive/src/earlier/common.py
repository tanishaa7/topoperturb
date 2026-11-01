"""Shared loading + PH helpers for the Phase 1 analyses."""
from pathlib import Path
import argparse
import os
import sys
import gzip
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from ripser import ripser
from persim import wasserstein

SEED = 42
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RESULTS = Path(os.environ.get("TOPOPERTURB_RESULTS", ROOT / "results/local"))
OUT = RESULTS / "phase1"
CACHE = DATA / "expression_k562.npz"
PILOT_CSV = DATA / "pilot.csv"
sys.path.insert(0, str(ROOT / "src"))
from prepare_data import load_cache


def configure_run(description, mapper=False):
    """Parse only output paths; retain each historical experiment's parameters."""
    global RESULTS, OUT
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--output-root", type=Path, default=RESULTS)
    if mapper:
        parser.add_argument("--skip-mapper", action="store_true", help="omit optional Mapper HTML views")
    args = parser.parse_args()
    RESULTS = args.output_root.resolve()
    OUT = RESULTS / "phase1"
    os.environ["TOPOPERTURB_RESULTS"] = str(RESULTS)
    for phase in ("pilot", "initial_ph", "phase1", "phase2"):
        (RESULTS / phase).mkdir(parents=True, exist_ok=True)
    return args


def load():
    """Controls + the eight pilot groups (never 'excluded'), cached as npz."""
    pilot = pd.read_csv(PILOT_CSV)
    if CACHE.exists():
        X, labels, genes = load_cache(CACHE, set(pilot.gene_id) | {"non-targeting"})
        return X, labels, genes, pilot
    from huggingface_hub import hf_hub_download
    path = hf_hub_download(repo_id="huggingbahl21/causalbench-k562",
                           filename="real/causalbench/data/causalbench_k562.tsv.gz",
                           repo_type="dataset")
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n\r").split("\t")
    assert header[-1] == "intervention"
    gene_cols = header[:-1]
    keep = set(pilot.gene_id) | {"non-targeting"}
    dtypes = {g: np.float32 for g in gene_cols}
    dtypes["intervention"] = "string"
    parts = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", dtype=dtypes, chunksize=5000):
        sub = chunk.loc[chunk["intervention"].isin(keep)]
        if len(sub):
            parts.append(sub)
    data = pd.concat(parts, ignore_index=True)
    X = data.loc[:, gene_cols].to_numpy(np.float32)
    lab = data["intervention"].astype(str).to_numpy()
    assert "excluded" not in set(lab)
    np.savez_compressed(CACHE, X=X, lab=lab.astype("U32"), genes=np.array(gene_cols, dtype="U32"))
    return X, lab, gene_cols, pilot


def project(X, n_pcs):
    """Fixed PCA basis fit once on all retained cells."""
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=SEED)
    return pca.fit_transform(X).astype(np.float32)


def dgms_of(points):
    d = ripser(points, maxdim=1)["dgms"]
    h0 = d[0][np.isfinite(d[0][:, 1])]           # drop the infinite H0 bar
    h1 = d[1] if len(d[1]) else np.empty((0, 2))
    return h0, h1


def wass(a, b):
    a = a if len(a) else np.empty((0, 2))
    b = b if len(b) else np.empty((0, 2))
    if len(a) == 0 and len(b) == 0:
        return 0.0
    return float(wasserstein(a, b))


def pair_distance(Z, idx_a, idx_b):
    """(H0, H1) Wasserstein between the diagrams of two point clouds."""
    a0, a1 = dgms_of(Z[idx_a])
    b0, b1 = dgms_of(Z[idx_b])
    return wass(a0, b0), wass(a1, b1)


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        prev = min(prev, p[i] * m / (m - rank))
        adj[i] = prev
    return np.clip(adj, 0, 1)
