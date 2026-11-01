"""Phase 1.3 — what do the H1 changes represent?

Per repeat (30 PCs, 150 cells): one control cloud, a second disjoint control
cloud (null) and one perturbed cloud. We compare
  raw        : H1 Wasserstein as before
  scaled     : each cloud divided by its own mean pairwise distance
               (removes global dispersion / scale)
  prominent  : only H1 features with persistence > T
  short      : only H1 features with persistence <= T
where T = 95th percentile of H1 persistence pooled over control clouds.
Representative cycles are extracted for the most persistent loop.
"""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import shortest_path
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from ripser import ripser
from common import load, project, wass, OUT, SEED

N, N_PCS, R = 150, 30, 30
X, lab, genes, pilot = load()
Z = project(X, N_PCS)
ctrl = np.where(lab == "non-targeting")[0]
sym = dict(zip(pilot.gene_id, pilot.gene_symbol))
by_sym = {v: k for k, v in sym.items()}
GENES = ["NCBP2", "POLR2B", "ZMAT2", "RPL21", "CPSF3", "GFM1", "DLD", "POLRMT"]
LEAD = ["NCBP2", "POLR2B", "ZMAT2"]


def h1(points):
    d = ripser(points, maxdim=1)["dgms"][1]
    return d if len(d) else np.empty((0, 2))


def pers(d):
    return d[:, 1] - d[:, 0]


def summary(points, d):
    p = pers(d)
    nn = NearestNeighbors(n_neighbors=6).fit(points)
    knn = nn.kneighbors(points)[0][:, 1:].mean(1)
    pr = p / p.sum() if p.sum() > 0 else p
    return {"n_H1": len(d), "total_pers": p.sum(), "max_pers": p.max() if len(p) else 0.0,
            "mean_birth": d[:, 0].mean() if len(d) else np.nan,
            "entropy": float(-(pr * np.log(pr + 1e-12)).sum()),
            "median_knn5": float(np.median(knn)), "mean_pairdist": float(pdist(points).mean())}


def control_threshold():
    rng = np.random.default_rng(SEED + 7)
    allp = []
    for _ in range(R):
        allp.append(pers(h1(Z[rng.choice(ctrl, N, replace=False)])))
    return float(np.quantile(np.concatenate(allp), 0.95))


T = None  # set in main and passed to workers


def run_gene(args):
    gsym, T = args
    pc = np.where(lab == by_sym[gsym])[0]
    rng = np.random.default_rng(SEED + 100 + GENES.index(gsym))
    rows = []
    for k in range(R):
        q = rng.choice(ctrl, 2 * N, replace=False)
        c, c2 = Z[q[:N]], Z[q[N:]]
        p = Z[rng.choice(pc, N, replace=False)]
        dc, dc2, dp = h1(c), h1(c2), h1(p)
        # scale-normalised clouds
        cs, c2s, ps = (a / pdist(a).mean() for a in (c, c2, p))
        dcs, dc2s, dps = h1(cs), h1(c2s), h1(ps)
        prom = lambda d: d[pers(d) > T]
        short = lambda d: d[pers(d) <= T]
        row = {"gene_symbol": gsym, "rep": k,
               "raw_obs": wass(dc, dp), "raw_null": wass(dc, dc2),
               "scaled_obs": wass(dcs, dps), "scaled_null": wass(dcs, dc2s),
               "prom_obs": wass(prom(dc), prom(dp)), "prom_null": wass(prom(dc), prom(dc2)),
               "short_obs": wass(short(dc), short(dp)), "short_null": wass(short(dc), short(dc2)),
               "n_prom_pert": int((pers(dp) > T).sum()), "n_prom_ctrl": int((pers(dc) > T).sum())}
        for pre, pts, d in (("pert_", p, dp), ("ctrl_", c, dc)):
            for key, v in summary(pts, d).items():
                row[pre + key] = v
        rows.append(row)
    return rows


def representative_cycle(points):
    """Most persistent H1 class -> an explicit loop of cells.
    Take the cocycle edge whose length matches the birth radius, then close it
    with the shortest path through edges shorter than the birth radius."""
    res = ripser(points, maxdim=1, do_cocycles=True)
    d = res["dgms"][1]
    if not len(d):
        return None
    i = int(np.argmax(pers(d)))
    b, de = d[i]
    D = squareform(pdist(points))
    edges = res["cocycles"][1][i][:, :2].astype(int)
    lens = D[edges[:, 0], edges[:, 1]]
    u, v = edges[np.argmin(np.abs(lens - b))]
    G = np.where(D < b - 1e-9, D, 0.0)
    G[u, v] = G[v, u] = 0.0
    dist, pred = shortest_path(G, directed=False, indices=u, return_predecessors=True)
    if not np.isfinite(dist[v]):
        return None
    path = [v]
    while path[-1] != u:
        path.append(pred[path[-1]])
    loop = path[::-1] + [u]   # u ... v, then the closing edge back to u
    return {"loop": loop, "birth": float(b), "death": float(de), "pers": float(de - b),
            "n_cells": len(path), "loop_len": float(sum(D[a, c] for a, c in zip(loop[:-1], loop[1:])))}


if __name__ == "__main__":
    t = time.time()
    T = control_threshold()
    print(f"Prominence threshold T (control H1 persistence p95) = {T:.3f}")
    with Pool(8) as pool:
        out = pool.map(run_gene, [(g, T) for g in GENES])
    df = pd.DataFrame([r for rows in out for r in rows])
    df.to_csv(OUT / "interpretation_repeats.csv", index=False)

    summ = []
    for g, s in df.groupby("gene_symbol", sort=False):
        r = {"gene_symbol": g}
        for m in ("raw", "scaled", "prom", "short"):
            r[f"{m}_ratio"] = s[f"{m}_obs"].mean() / s[f"{m}_null"].mean()
            r[f"{m}_win_rate"] = float((s[f"{m}_obs"] > s[f"{m}_null"]).mean())
        for key in ("n_H1", "total_pers", "max_pers", "entropy", "median_knn5", "mean_pairdist", "mean_birth"):
            r[f"{key}_pert"] = s[f"pert_{key}"].mean()
            r[f"{key}_ctrl"] = s[f"ctrl_{key}"].mean()
        r["n_prom_pert"] = s.n_prom_pert.mean(); r["n_prom_ctrl"] = s.n_prom_ctrl.mean()
        summ.append(r)
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "interpretation_summary.csv", index=False)

    # ---- persistence diagrams: control vs lead genes, H1, greyscale ----
    rng = np.random.default_rng(SEED + 55)
    c_idx = rng.choice(ctrl, N, replace=False)
    dc = h1(Z[c_idx])
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True, sharey=True)
    lim, lo = 0, np.inf
    clouds = {}
    for ax, g in zip(axes, LEAD):
        p_idx = rng.choice(np.where(lab == by_sym[g])[0], N, replace=False)
        clouds[g] = p_idx
        dp = h1(Z[p_idx])
        lim = max(lim, dc.max(), dp.max()); lo = min(lo, dc[:, 0].min(), dp[:, 0].min())
        ax.scatter(dc[:, 0], dc[:, 1], s=18, facecolors="none", edgecolors="0.45", label="control")
        ax.scatter(dp[:, 0], dp[:, 1], s=16, color="black", marker="x", label=g)
        ax.set_title(f"H1: control vs {g}")
        ax.set_xlabel("birth")
    for ax in axes:
        xs = np.array([lo * 0.95, lim * 1.05])
        ax.plot(xs, xs, color="black", lw=0.8)
        ax.plot(xs, xs + T, color="black", lw=0.8, ls="--", label=f"prominence T = {T:.2f}")
        ax.set_xlim(lo * 0.95, lim * 1.05); ax.set_ylim(lo * 0.95, lim * 1.05)
        ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    axes[0].set_ylabel("death")
    fig.tight_layout()
    fig.savefig(OUT / "h1_diagrams_lead_genes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- representative cycles (most persistent loop), greyscale ----
    cyc_rows = []
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, (name, idx) in zip(axes, [("control", c_idx)] + [(g, clouds[g]) for g in LEAD]):
        pts = Z[idx]
        cyc = representative_cycle(pts)
        if cyc is None:
            ax.set_title(f"{name}: no loop"); continue
        loop = cyc["loop"]
        local = PCA(n_components=2, random_state=SEED).fit(pts[loop[:-1]])
        P = local.transform(pts)
        ax.scatter(P[:, 0], P[:, 1], s=6, color="0.8")
        ax.plot(P[loop, 0], P[loop, 1], color="black", lw=1.2)
        ax.scatter(P[loop, 0], P[loop, 1], s=18, color="black")
        ax.set_title(f"{name}\nloop of {cyc['n_cells']} cells, persistence {cyc['pers']:.2f}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        cyc_rows.append({"group": name, **{k: v for k, v in cyc.items() if k != "loop"}})
    fig.suptitle("Most persistent H1 loop in one 150-cell sample (cells shown in the loop's own 2-D plane)", y=1.03)
    fig.tight_layout()
    fig.savefig(OUT / "representative_cycles.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(cyc_rows).to_csv(OUT / "representative_cycles.csv", index=False)

    pd.set_option("display.width", 200)
    print(summ[["gene_symbol", "raw_ratio", "scaled_ratio", "prom_ratio", "short_ratio",
                "raw_win_rate", "scaled_win_rate", "prom_win_rate"]].round(2).to_string(index=False))
    print()
    print(summ[["gene_symbol", "n_H1_pert", "n_H1_ctrl", "n_prom_pert", "n_prom_ctrl",
                "max_pers_pert", "max_pers_ctrl", "mean_pairdist_pert", "mean_pairdist_ctrl",
                "median_knn5_pert", "median_knn5_ctrl"]].round(2).to_string(index=False))
    print()
    print(pd.DataFrame(cyc_rows).round(3).to_string(index=False))
    print(f"\n{time.time()-t:.0f}s")
