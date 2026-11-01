"""Phase 2.4 — independent RPE1 validation.

RPE1 is analysed on its own: own PCA (30 PCs, fit on RPE1 controls + tested
genes), own controls, same statistics as K562.
Shared genes: ZMAT2, CPSF3, POLRMT, GFM1 (POLR2B, RPL21, NCBP2, DLD are not
perturbed in RPE1). Positive control: the RPE1 knockdown (>=150 cells, not
one of the four) with the largest centroid shift in full gene space —
chosen for its large shift, so it is a pipeline check, not a test.
Per gene, one 150-vs-150 draw and 2,000 label permutations of:
  raw H1 Wasserstein, scale-normalised H1 Wasserstein, energy distance.
Normalising a cloud by its mean pairwise distance scales its Rips diagram
by exactly the same factor, so both H1 statistics share one ripser call.
"""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time, gzip
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from ripser import ripser
from common import project, wass, bh, SEED, DATA, RESULTS, load_cache

N, N_PCS, N_PERM, N_CNULL = 150, 30, 2000, 500
SHARED = {"ENSG00000146007": "ZMAT2", "ENSG00000119203": "CPSF3",
          "ENSG00000099821": "POLRMT", "ENSG00000168827": "GFM1"}
CACHE = DATA / "expression_rpe1.npz"
OUT2 = RESULTS / "phase2"
OUT2.mkdir(exist_ok=True)


def load_rpe1():
    if CACHE.exists():
        return load_cache(CACHE)
    from huggingface_hub import hf_hub_download
    path = hf_hub_download(repo_id="huggingbahl21/causalbench-rpe1",
                           filename="real/causalbench/data/causalbench_rpe1.tsv.gz", repo_type="dataset")
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n\r").split("\t")
    assert header[-1] == "intervention"
    gene_cols = header[:-1]
    dtypes = {g: np.float32 for g in gene_cols}
    dtypes["intervention"] = "string"
    parts = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", dtype=dtypes, chunksize=10000):
        parts.append(chunk.loc[chunk["intervention"] != "excluded"])
    data = pd.concat(parts, ignore_index=True)
    X = data.loc[:, gene_cols].to_numpy(np.float32)
    lab = data["intervention"].astype(str).to_numpy()
    np.savez_compressed(CACHE, X=X, lab=lab.astype("U32"), genes=np.array(gene_cols, dtype="U32"))
    return X, lab, gene_cols


Xfull, labfull, gene_cols = load_rpe1()
ctrl_mask = labfull == "non-targeting"
counts = pd.Series(labfull[~ctrl_mask]).value_counts()
mu_c = Xfull[ctrl_mask].mean(0)
screen = pd.DataFrame([{"gene_id": g, "n_cells": int(n),
                        "centroid_shift_fullgene": float(np.linalg.norm(Xfull[labfull == g].mean(0) - mu_c))}
                       for g, n in counts.items() if n >= N and g not in SHARED])
screen = screen.sort_values("centroid_shift_fullgene", ascending=False)
POS = screen.iloc[0].gene_id
TEST = {**SHARED, POS: "positive_control"}

keep = np.isin(labfull, list(TEST) + ["non-targeting"])
X, lab = Xfull[keep], labfull[keep]
Z = project(X, N_PCS)
ctrl = np.where(lab == "non-targeting")[0]


def h1_and_scale(P):
    d = ripser(P, maxdim=1)["dgms"][1]
    return (d if len(d) else np.empty((0, 2))), pdist(P).mean()


def energy(A, B):
    return 2 * np.mean(np.linalg.norm(A[:, None] - B[None], axis=2)) - pdist(A).mean() - pdist(B).mean()


def stat3(A, B):
    (da, sa), (db, sb) = h1_and_scale(A), h1_and_scale(B)
    return wass(da, db), wass(da / sa, db / sb), energy(A, B), sb / sa


def run_gene(k):
    g = list(TEST)[k]
    rng = np.random.default_rng(SEED + 2000 * (k + 1))
    pc = np.where(lab == g)[0]
    c = rng.choice(ctrl, N, replace=False)
    p = rng.choice(pc, N, replace=False) if len(pc) > N else pc
    raw, nrm, en, dil = stat3(Z[c], Z[p])
    cs = np.linalg.norm(Z[c].mean(0) - Z[p].mean(0))
    cnull = np.empty(N_CNULL)
    for j in range(N_CNULL):
        q = rng.choice(ctrl, 2 * N, replace=False)
        cnull[j] = np.linalg.norm(Z[q[:N]].mean(0) - Z[q[N:]].mean(0))
    pool = np.concatenate([c, p])
    perm = np.empty((N_PERM, 3))
    for j in range(N_PERM):
        q = rng.permutation(pool)
        perm[j] = stat3(Z[q[:N]], Z[q[N:]])[:3]
    pv = lambda o, col: (1 + (perm[:, col] >= o).sum()) / (1 + N_PERM)
    return {"gene_id": g, "label": TEST[g], "n_cells_rpe1": len(pc),
            "centroid_shift": cs, "centroid_null_p95": np.quantile(cnull, 0.95),
            "p_centroid": (1 + (cnull >= cs).sum()) / (1 + N_CNULL),
            "cloud_size_ratio": dil,
            "raw_H1": raw, "raw_H1_vs_perm": raw / perm[:, 0].mean(), "p_raw_H1": pv(raw, 0),
            "norm_H1": nrm, "norm_H1_vs_perm": nrm / perm[:, 1].mean(), "p_norm_H1": pv(nrm, 1),
            "energy": en, "energy_vs_perm": en / perm[:, 2].mean(), "p_energy": pv(en, 2)}


if __name__ == "__main__":
    t = time.time()
    print(f"RPE1 retained (non-excluded): {Xfull.shape}; controls {int(ctrl_mask.sum())}")
    print("Positive-control screen (top 5 by full-gene centroid shift, >=150 cells, excluding the 4 shared):")
    print(screen.head(5).round(3).to_string(index=False))
    screen.to_csv(OUT2 / "rpe1_centroid_screen.csv", index=False)
    with Pool(5) as pool:
        rows = pool.map(run_gene, range(len(TEST)))
    res = pd.DataFrame(rows)
    sh = res.label != "positive_control"
    for col in ("p_centroid", "p_raw_H1", "p_norm_H1", "p_energy"):
        res.loc[sh, col.replace("p_", "q_")] = bh(res.loc[sh, col].values)
    res.to_csv(OUT2 / "rpe1_validation.csv", index=False)
    pd.set_option("display.width", 250)
    print(res.round(4).to_string(index=False))
    print(f"\n{time.time()-t:.0f}s")
