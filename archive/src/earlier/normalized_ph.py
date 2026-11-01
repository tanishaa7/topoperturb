"""Phase 2.1 + 2.3 — scale-normalised H1 permutation tests (K562).

Every cloud (observed and each permuted half) is turned into a distance
matrix, divided by its own mean off-diagonal distance, then given to ripser.
The statistic is the H1 Wasserstein distance between the two diagrams.

Metrics
  euclid    : Euclidean distance in PCA space
  effres    : corrected effective resistance on a symmetric kNN graph (k=15)
  diffusion : diffusion distance on the same graph (t=8)
(spectral distances after Damrich, Berens & Kobak, NeurIPS 2024)

Tasks
  primary  : euclid, 30 PCs, 150 cells, all 8 genes, 2,000 permutations
  seed     : euclid, 30 PCs, 150 cells, focus genes, 2 new draws, 1,000 permutations
  grid     : euclid, focus genes, PCs {10,20,30,50} x cells {100,150,200}, 200 permutations
  spectral : effres + diffusion, 30 PCs, 150 cells, all 8 genes, 500 permutations
"""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time, csv
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from ripser import ripser
from common import load, project, wass, bh, OUT, SEED

FOCUS = ["POLR2B", "RPL21", "NCBP2"]
DIMS, NS = [10, 20, 30, 50], [100, 150, 200]

X, lab, genes, pilot = load()
ctrl = np.where(lab == "non-targeting")[0]
ids = pilot.gene_id.tolist()
sym = dict(zip(pilot.gene_id, pilot.gene_symbol))
by_sym = {v: k for k, v in sym.items()}
_Z = {}


def Zd(d):
    if d not in _Z:
        _Z[d] = project(X, d)
    return _Z[d]


from methods import dist_effres, dist_diffusion


def dist_euclid(P):
    return squareform(pdist(P))


METRICS = {"euclid": dist_euclid, "effres": dist_effres, "diffusion": dist_diffusion}


def h1_norm(D):
    D = D / D[np.triu_indices(len(D), 1)].mean()
    d = ripser(D, maxdim=1, distance_matrix=True)["dgms"][1]
    return d if len(d) else np.empty((0, 2))


def run(task):
    kind, gsym, d, n, seed_k, metric, nperm = task
    t0 = time.time()
    Z = Zd(d)
    g = by_sym[gsym]
    i = ids.index(g)
    # seed_k = 0 with 30 PCs / 150 cells reproduces the Phase 1.1 draw
    rng = np.random.default_rng(SEED + 1000 * (i + 1) + 100000 * seed_k)
    c = rng.choice(ctrl, n, replace=False)
    pc = np.where(lab == g)[0]
    p = rng.choice(pc, n, replace=False) if len(pc) > n else pc
    f = METRICS[metric]
    obs = wass(h1_norm(f(Z[c])), h1_norm(f(Z[p])))
    pool = np.concatenate([c, p])
    perm = np.empty(nperm)
    for k in range(nperm):
        q = rng.permutation(pool)
        perm[k] = wass(h1_norm(f(Z[q[:n]])), h1_norm(f(Z[q[n:]])))
    n_ge = int((perm >= obs).sum())
    return {"kind": kind, "gene_symbol": gsym, "n_pcs": d, "n_cells": n, "seed": seed_k,
            "metric": metric, "n_perm": nperm, "obs": obs, "perm_mean": perm.mean(),
            "perm_p95": np.quantile(perm, 0.95), "effect_vs_perm": obs / perm.mean(),
            "n_ge": n_ge, "p": (1 + n_ge) / (1 + nperm),
            "n_available": len(pc), "perturbed_side_fixed": len(pc) <= n,
            "seconds": round(time.time() - t0)}


if __name__ == "__main__":
    allsym = [sym[g] for g in ids]
    tasks = [("primary", g, 30, 150, 0, "euclid", 2000) for g in allsym]
    tasks += [("seed", g, 30, 150, k, "euclid", 1000) for g in FOCUS for k in (1, 2)]
    tasks += [("spectral", g, 30, 150, 0, m, 500) for g in allsym for m in ("effres", "diffusion")]
    tasks += [("grid", g, d, n, 0, "euclid", 200) for g in FOCUS for d in DIMS for n in NS]
    tasks.sort(key=lambda t: -t[6] * t[3] ** 2)
    out = OUT.parent / "phase2"
    out.mkdir(exist_ok=True)
    path = out / "p2_ph_raw_results.csv"
    t = time.time()
    with Pool(8) as pool, open(path, "w", newline="") as fh:
        w = None
        for r in pool.imap_unordered(run, tasks, chunksize=1):
            if w is None:
                w = csv.DictWriter(fh, fieldnames=list(r)); w.writeheader()
            w.writerow(r); fh.flush()
            print(f"[{time.time()-t:5.0f}s] {r['kind']:8s} {r['gene_symbol']:7s} {r['metric']:9s} "
                  f"{r['n_pcs']}PC/{r['n_cells']} seed{r['seed']}  obs={r['obs']:.4f} "
                  f"perm_mean={r['perm_mean']:.4f}  p={r['p']:.4f}", flush=True)
    print(f"done {time.time()-t:.0f}s")
