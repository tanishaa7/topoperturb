"""Phase 1.1 — 2,000-permutation test for all eight genes (30 PCs, 150 cells/group)."""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from common import load, project, pair_distance, bh, OUT, SEED

N, N_PCS, N_PERM = 150, 30, 2000

X, lab, genes, pilot = load()
Z = project(X, N_PCS)
ctrl = np.where(lab == "non-targeting")[0]
ids = pilot.gene_id.tolist()
sym = dict(zip(pilot.gene_id, pilot.gene_symbol))


def run_gene(i):
    g = ids[i]
    rng = np.random.default_rng(SEED + 1000 * (i + 1))   # fixed per-gene seed
    c = rng.choice(ctrl, N, replace=False)
    p = rng.choice(np.where(lab == g)[0], N, replace=False)
    obs0, obs1 = pair_distance(Z, c, p)
    pool = np.concatenate([c, p])
    perm0 = np.empty(N_PERM); perm1 = np.empty(N_PERM)
    for k in range(N_PERM):
        q = rng.permutation(pool)
        perm0[k], perm1[k] = pair_distance(Z, q[:N], q[N:])
    return {"gene_id": g, "gene_symbol": sym[g],
            "wass_H0": obs0, "wass_H1": obs1,
            "perm_H1_mean": perm1.mean(), "perm_H1_p95": np.quantile(perm1, 0.95),
            "n_perm_H1_ge_obs": int((perm1 >= obs1).sum()),
            "perm_p_H0": (1 + (perm0 >= obs0).sum()) / (1 + N_PERM),
            "perm_p_H1": (1 + (perm1 >= obs1).sum()) / (1 + N_PERM)}


if __name__ == "__main__":
    t = time.time()
    with Pool(8) as pool:
        rows = pool.map(run_gene, range(len(ids)))
    res = pd.DataFrame(rows)
    res["bh_q_H0"] = bh(res.perm_p_H0)
    res["bh_q_H1"] = bh(res.perm_p_H1)
    res["significant_H1_q05"] = res.bh_q_H1 < 0.05
    res = res.sort_values("perm_p_H1")
    res.to_csv(OUT / "perm2000_significance.csv", index=False)
    print(res[["gene_symbol", "wass_H1", "perm_H1_mean", "perm_H1_p95", "n_perm_H1_ge_obs",
               "perm_p_H1", "bh_q_H1", "perm_p_H0", "bh_q_H0"]].round(4).to_string(index=False))
    print(f"\n{time.time()-t:.0f}s")
