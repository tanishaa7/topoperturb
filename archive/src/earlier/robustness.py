"""Phase 1.2 — H1 robustness across PCA dimension x cells per group."""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time, itertools
from multiprocessing import Pool
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load, project, pair_distance, OUT, SEED

DIMS = [10, 20, 30, 50]
NS = [100, 150, 200]
N_REPEAT = 30
N_NULL = 200
N_BOOT = 2000

X, lab, genes, pilot = load()
ctrl = np.where(lab == "non-targeting")[0]
ids = pilot.gene_id.tolist()
sym = dict(zip(pilot.gene_id, pilot.gene_symbol))
LEAD = ["NCBP2", "POLR2B", "ZMAT2"]


def run_config(cfg):
    d, n = cfg
    Z = project(X, d)
    rng = np.random.default_rng(SEED + 10 * d + n)
    null = np.empty(N_NULL)
    for k in range(N_NULL):
        q = rng.choice(ctrl, 2 * n, replace=False)
        null[k] = pair_distance(Z, q[:n], q[n:])[1]
    rows = []
    for g in ids:
        pc = np.where(lab == g)[0]
        rep = np.empty(N_REPEAT)
        for k in range(N_REPEAT):
            c = rng.choice(ctrl, n, replace=False)
            p = rng.choice(pc, n, replace=False) if len(pc) > n else pc
            rep[k] = pair_distance(Z, c, p)[1]
        boot = rng.choice(rep, (N_BOOT, N_REPEAT), replace=True).mean(1)
        nm = null.mean()
        rows.append({"n_pcs": d, "n_cells": n, "gene_id": g, "gene_symbol": sym[g],
                     "n_available": len(pc), "perturbed_side_fixed": len(pc) <= n,
                     "H1_mean": rep.mean(), "H1_sd": rep.std(ddof=1),
                     "H1_mean_ci_lo": np.quantile(boot, 0.025), "H1_mean_ci_hi": np.quantile(boot, 0.975),
                     "null_mean": nm, "null_sd": null.std(ddof=1), "null_p95": np.quantile(null, 0.95),
                     "ratio_to_null": rep.mean() / nm,
                     "ratio_ci_lo": np.quantile(boot, 0.025) / nm, "ratio_ci_hi": np.quantile(boot, 0.975) / nm,
                     "z_vs_null": (rep.mean() - nm) / null.std(ddof=1),
                     "detect_rate": float((rep > np.quantile(null, 0.95)).mean())})
    return rows


if __name__ == "__main__":
    t = time.time()
    with Pool(8) as pool:
        out = pool.map(run_config, list(itertools.product(DIMS, NS)))
    res = pd.DataFrame([r for rows in out for r in rows])
    res.to_csv(OUT / "robustness_grid.csv", index=False)

    # ---- heatmap: ratio to matched control null (greyscale) ----
    res["config"] = res.n_pcs.astype(str) + " PCs / " + res.n_cells.astype(str)
    order_cfg = [f"{d} PCs / {n}" for d in DIMS for n in NS]
    piv = res.pivot(index="gene_symbol", columns="config", values="ratio_to_null")[order_cfg]
    det = res.pivot(index="gene_symbol", columns="config", values="detect_rate")[order_cfg]
    piv = piv.loc[piv.mean(1).sort_values(ascending=False).index]
    det = det.loc[piv.index]
    fig, ax = plt.subplots(figsize=(13, 4.6))
    im = ax.imshow(piv.values, cmap="Greys", vmin=0.6, vmax=3.2, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            ax.text(j, i, f"{v:.2f}" + ("*" if det.values[i, j] >= 0.8 else ""),
                    ha="center", va="center", fontsize=8,
                    color="white" if v > 2.0 else "black")
    ax.set_xticks(range(len(order_cfg)))
    ax.set_xticklabels([c.replace(" / ", "\n") + " cells" for c in order_cfg], fontsize=7.5)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
    for j in (3, 6, 9):
        ax.axvline(j - 0.5, color="black", lw=1)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("H1 Wasserstein / control-null mean")
    ax.set_title("H1 topological effect relative to matched control null "
                 "(* = above null 95th pct in ≥80% of 30 repeats)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "robustness_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- CI plot: ratio + bootstrap 95% CI, one panel per cell count ----
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)
    markers = {"NCBP2": "o", "POLR2B": "s", "ZMAT2": "^"}
    for ax, n in zip(axes, NS):
        sub = res[res.n_cells == n]
        for g, gs in sub.groupby("gene_symbol"):
            gs = gs.sort_values("n_pcs")
            lead = g in LEAD
            off = {"NCBP2": -1.2, "POLR2B": 0, "ZMAT2": 1.2}.get(g, 0)
            ax.errorbar(gs.n_pcs + off, gs.ratio_to_null,
                        yerr=[gs.ratio_to_null - gs.ratio_ci_lo, gs.ratio_ci_hi - gs.ratio_to_null],
                        color="black" if lead else "0.75", marker=markers.get(g, "."),
                        ms=5 if lead else 3, lw=1.3 if lead else 0.8, capsize=2,
                        label=g if lead else None, zorder=3 if lead else 1)
        ax.axhline(1.0, color="black", ls=":", lw=1)
        ax.set_title(f"{n} cells per group"); ax.set_xlabel("PCA components")
        ax.set_xticks(DIMS)
    axes[0].set_ylabel("H1 effect / control-null mean (95% CI)")
    axes[0].plot([], [], color="0.75", label="other 5 genes"); axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("Robustness of the H1 effect across PCA dimension and sample size", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "robustness_ci.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    show = res[res.gene_symbol.isin(LEAD + ["RPL21", "CPSF3", "GFM1", "DLD", "POLRMT"])]
    print(show[["gene_symbol", "n_pcs", "n_cells", "ratio_to_null", "ratio_ci_lo", "ratio_ci_hi",
                "detect_rate", "z_vs_null"]].round(2).to_string(index=False))
    print(f"\n{time.time()-t:.0f}s")
