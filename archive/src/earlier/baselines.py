"""Phase 2.2 — conventional statistics vs persistent homology (K562).

Same draw as Phase 1.1 / Phase 2 primary (30 PCs, 150 cells/group) and
2,000 label permutations. Each half is transformed separately, exactly as
the PH statistic is:
  raw      : as is
  centered : each cloud minus its own mean
  cs       : centered, then divided by its own mean pairwise distance
Measures (two-sided where it is a difference of summaries):
  centroid shift, dispersion difference, energy distance, MMD (RBF,
  bandwidth = median pooled distance of the observed draw), covariance
  Frobenius difference, participation-ratio difference (anisotropy).
Gaussian surrogates: clouds drawn from N(mean, cov) of each full cell pool,
scale-normalised, H1 Wasserstein vs a surrogate control-null (30 draws).
"""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from ripser import ripser
from common import load, project, wass, bh, OUT, SEED

N, N_PCS, N_PERM, R_SURR = 150, 30, 2000, 30
X, lab, genes, pilot = load()
Z = project(X, N_PCS)
ctrl = np.where(lab == "non-targeting")[0]
ids = pilot.gene_id.tolist()
sym = dict(zip(pilot.gene_id, pilot.gene_symbol))
OUT2 = OUT.parent / "phase2"
OUT2.mkdir(exist_ok=True)


def transform(P, mode):
    if mode == "raw":
        return P
    P = P - P.mean(0)
    if mode == "cs":
        P = P / pdist(P).mean()
    return P


def energy(A, B):
    return 2 * np.mean(np.linalg.norm(A[:, None] - B[None], axis=2)) - pdist(A).mean() - pdist(B).mean()


def mmd2(A, B, sigma):
    k = lambda D: np.exp(-D ** 2 / (2 * sigma ** 2))
    Kaa = k(squareform(pdist(A))); Kbb = k(squareform(pdist(B)))
    Kab = k(np.linalg.norm(A[:, None] - B[None], axis=2))
    n = len(A)
    return (Kaa.sum() - n) / (n * (n - 1)) + (Kbb.sum() - n) / (n * (n - 1)) - 2 * Kab.mean()


def pr(C):
    return np.trace(C) ** 2 / np.trace(C @ C)


def stats(A, B, sig):
    out = {}
    out["raw:centroid"] = np.linalg.norm(A.mean(0) - B.mean(0))
    rad = lambda P: np.linalg.norm(P - P.mean(0), axis=1).mean()
    out["raw:dispersion"] = abs(rad(A) - rad(B))
    for mode in ("raw", "centered", "cs"):
        a, b = transform(A, mode), transform(B, mode)
        out[f"{mode}:energy"] = energy(a, b)
        out[f"{mode}:mmd"] = mmd2(a, b, sig[mode])
        if mode != "centered":           # covariance does not change with centering
            out[f"{mode}:cov_frob"] = np.linalg.norm(np.cov(a.T) - np.cov(b.T))
    out["cs:anisotropy_PR"] = abs(pr(np.cov(A.T)) - pr(np.cov(B.T)))   # scale + shift invariant
    return out


def run_gene(i):
    g = ids[i]
    rng = np.random.default_rng(SEED + 1000 * (i + 1))      # same draw as Phase 1.1
    c = rng.choice(ctrl, N, replace=False)
    pc = np.where(lab == g)[0]
    p = rng.choice(pc, N, replace=False)
    A, B = Z[c], Z[p]
    sig = {m: np.median(pdist(np.vstack([transform(A, m), transform(B, m)]))) for m in ("raw", "centered", "cs")}
    obs = stats(A, B, sig)
    pool = np.concatenate([c, p])
    ge = {k: 0 for k in obs}
    pm = {k: 0.0 for k in obs}
    prng = np.random.default_rng(SEED + 555 + i)
    for _ in range(N_PERM):
        q = prng.permutation(pool)
        s = stats(Z[q[:N]], Z[q[N:]], sig)
        for k, v in s.items():
            ge[k] += v >= obs[k]
            pm[k] += v / N_PERM
    rows = [{"gene_symbol": sym[g], "measure": k, "obs": obs[k], "perm_mean": pm[k],
             "effect_vs_perm": obs[k] / pm[k] if pm[k] else np.nan,
             "p": (1 + ge[k]) / (1 + N_PERM)} for k in obs]

    # Gaussian surrogates with each pool's exact mean + covariance
    srng = np.random.default_rng(SEED + 777 + i)
    mc, Cc = Z[ctrl].mean(0), np.cov(Z[ctrl].T)
    mp, Cp = Z[pc].mean(0), np.cov(Z[pc].T)

    def h1n(P):
        P = P / pdist(P).mean()
        d = ripser(P, maxdim=1)["dgms"][1]
        return d if len(d) else np.empty((0, 2))
    so, sn = [], []
    for _ in range(R_SURR):
        gc = h1n(srng.multivariate_normal(mc, Cc, N))
        so.append(wass(gc, h1n(srng.multivariate_normal(mp, Cp, N))))
        sn.append(wass(gc, h1n(srng.multivariate_normal(mc, Cc, N))))
    so, sn = np.array(so), np.array(sn)
    surr = {"gene_symbol": sym[g], "surr_ratio": so.mean() / sn.mean(),
            "surr_wins": int((so > sn).sum()), "surr_obs_mean": so.mean(), "surr_null_mean": sn.mean()}
    return rows, surr


if __name__ == "__main__":
    t = time.time()
    with Pool(4) as pool:
        res = pool.map(run_gene, range(len(ids)))
    df = pd.DataFrame([r for rows, _ in res for r in rows])
    df["q"] = df.groupby("measure")["p"].transform(lambda s: bh(s.values))
    df.to_csv(OUT2 / "conventional_measures.csv", index=False)
    surr = pd.DataFrame([s for _, s in res])
    surr.to_csv(OUT2 / "gaussian_surrogates.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.pivot(index="gene_symbol", columns="measure", values="q").round(4).to_string())
    print()
    print(df.pivot(index="gene_symbol", columns="measure", values="effect_vs_perm").round(2).to_string())
    print()
    print(surr.round(3).to_string(index=False))
    print(f"\n{time.time()-t:.0f}s")
