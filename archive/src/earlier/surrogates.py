"""Phase 2.5 — Gaussian surrogate check for every distance metric (K562).

Per gene and metric, 30 draws (30 PCs, 150 cells) of
  real      : control c, disjoint control c2, perturbed p (real cells)
  surrogate : the same three clouds drawn from N(mean, cov) of each full pool
Effect per draw = W(c, p) - W(c, c2), H1 Wasserstein on mean-normalised
distance matrices (same normalisation as Phase 2.1 / 2.3).
If surrogates reproduce the real effect, the PH signal is covariance-level.
Real vs surrogate effects are compared with a one-sided Mann-Whitney test
(real > surrogate), BH across genes per metric.
"""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

import sys, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from common import load, project, wass, bh, SEED, RESULTS
from normalized_ph import METRICS, h1_norm

N, N_PCS, R = 150, 30, 30
X, lab, genes, pilot = load()
Z = project(X, N_PCS)
ctrl = np.where(lab == "non-targeting")[0]
ids = pilot.gene_id.tolist()
sym = dict(zip(pilot.gene_id, pilot.gene_symbol))
OUT2 = str(RESULTS / "phase2") + "/"


def run(args):
    i, metric = args
    g = ids[i]
    f = METRICS[metric]
    H = lambda P: h1_norm(f(P))
    rng = np.random.default_rng(SEED + 4000 + 10 * i + list(METRICS).index(metric))
    pc = np.where(lab == g)[0]
    mc, Cc = Z[ctrl].mean(0), np.cov(Z[ctrl].T)
    mp, Cp = Z[pc].mean(0), np.cov(Z[pc].T)
    rows = []
    for k in range(R):
        q = rng.choice(ctrl, 2 * N, replace=False)
        p = rng.choice(pc, N, replace=False)
        hc = H(Z[q[:N]])
        ro, rn = wass(hc, H(Z[p])), wass(hc, H(Z[q[N:]]))
        sc = H(rng.multivariate_normal(mc, Cc, N))
        so = wass(sc, H(rng.multivariate_normal(mp, Cp, N)))
        sn = wass(sc, H(rng.multivariate_normal(mc, Cc, N)))
        rows.append({"gene_symbol": sym[g], "metric": metric, "draw": k,
                     "real_obs": ro, "real_null": rn, "surr_obs": so, "surr_null": sn})
    return rows


if __name__ == "__main__":
    t = time.time()
    with Pool(8) as pool:
        out = pool.map(run, [(i, m) for i in range(len(ids)) for m in METRICS])
    df = pd.DataFrame([r for rows in out for r in rows])
    df.to_csv(OUT2 + "surrogate_all_metrics_draws.csv", index=False)
    summ = []
    for (g, m), s in df.groupby(["gene_symbol", "metric"]):
        re, se = s.real_obs - s.real_null, s.surr_obs - s.surr_null
        summ.append({"gene_symbol": g, "metric": m,
                     "real_ratio": s.real_obs.mean() / s.real_null.mean(),
                     "real_wins": int((re > 0).sum()),
                     "surr_ratio": s.surr_obs.mean() / s.surr_null.mean(),
                     "surr_wins": int((se > 0).sum()),
                     "real_minus_surr_effect": re.mean() - se.mean(),
                     "p_real_gt_surr": mannwhitneyu(re, se, alternative="greater").pvalue})
    summ = pd.DataFrame(summ)
    summ["q_real_gt_surr"] = summ.groupby("metric")["p_real_gt_surr"].transform(lambda x: bh(x.values))
    summ.to_csv(OUT2 + "surrogate_all_metrics_summary.csv", index=False)
    pd.set_option("display.width", 200)
    for m in METRICS:
        print(f"== {m} ==")
        print(summ[summ.metric == m].sort_values("real_ratio", ascending=False)
              .drop(columns="metric").round(3).to_string(index=False))
    print(f"\n{time.time()-t:.0f}s")
