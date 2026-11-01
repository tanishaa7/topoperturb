"""Historical initial PH and target-gene omission checks.

Requires pilot.py outputs under the same --output-root. This original analysis
uses Python hash(gene) to seed draws. Set PYTHONHASHSEED for new repeatable runs;
the hash seed of the archived results was not recorded, so exact replay is not
promised. Main Phase 1/2 analyses use explicit numeric seeds instead.
"""

if __name__ == "__main__":
    from common import configure_run
    configure_run(__doc__)

from pathlib import Path
from common import DATA, RESULTS, load
import gzip, json, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from ripser import ripser
from persim import wasserstein, plot_diagrams

SEED = 42
N = 150                 # cells per group
N_PCS = 30
N_REPEAT = 20           # repeated-sampling draws for effect-size variation
N_PERM = 100            # permutations (feasibility, as requested)
N_CTRL_NULL = 200       # matched control-vs-control topological null draws
N_CENTROID_NULL = 500   # 150v150 centroid null draws
rng = np.random.default_rng(SEED)

OUT = RESULTS / "initial_ph"
OUT.mkdir(exist_ok=True)
DIAG_DIR = OUT / "diagrams"; DIAG_DIR.mkdir(exist_ok=True)

REPOS = {"K562": ("huggingbahl21/causalbench-k562", "real/causalbench/data/causalbench_k562.tsv.gz")}
PILOT = pd.read_csv(DATA / "pilot.csv")
SYM = dict(zip(PILOT.gene_id, PILOT.gene_symbol))
pilot_ids = PILOT.gene_id.tolist()
print("Pilot genes:", [(g, SYM[g]) for g in pilot_ids])
if not (RESULTS / "pilot/pilot_pca_effect_summary.csv").exists():
    raise FileNotFoundError("Run pilot.py first with the same --output-root")

# Reuse the recovered expression cache, preserving its row and gene order.
Xall, lab, gene_cols, _ = load()
gidx = {g: i for i, g in enumerate(gene_cols)}
print("Loaded:", Xall.shape, "controls:", int((lab == "non-targeting").sum()))

# ---- fixed PCA basis on all retained cells ----
pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=SEED)
Zall = pca.fit_transform(Xall).astype(np.float32)
ctrl_idx = np.where(lab == "non-targeting")[0]
pert_idx = {g: np.where(lab == g)[0] for g in pilot_ids}

def dgms_of(points):
    d = ripser(points, maxdim=1)["dgms"]
    h0 = d[0]; h0 = h0[np.isfinite(h0[:, 1])]          # drop the infinite H0 generator
    h1 = d[1] if len(d[1]) else np.empty((0, 2))
    return h0, h1

def wass(a, b):
    a = a if len(a) else np.empty((0, 2)); b = b if len(b) else np.empty((0, 2))
    if len(a) == 0 and len(b) == 0:
        return 0.0
    return float(wasserstein(a, b))

def pair_distance(Z, idx_a, idx_b):
    """H0/H1 Wasserstein between the diagrams of two point clouds."""
    a0, a1 = dgms_of(Z[idx_a]); b0, b1 = dgms_of(Z[idx_b])
    return wass(a0, b0), wass(a1, b1)

# ================= 1. corrected 150v150 CENTROID null =================
cshift_null = []
for _ in range(N_CENTROID_NULL):
    p = rng.permutation(ctrl_idx); a, b = p[:N], p[N:2*N]
    cshift_null.append(np.linalg.norm(Zall[a].mean(0) - Zall[b].mean(0)))
cshift_null = np.array(cshift_null)
cs_mean, cs_p95 = float(cshift_null.mean()), float(np.quantile(cshift_null, 0.95))
print(f"\n[150v150] centroid null: mean={cs_mean:.3f}, p95={cs_p95:.3f}")

# ================= 2. matched-control TOPOLOGICAL null (shared) =================
t0 = time.time()
null_h0, null_h1 = [], []
for _ in range(N_CTRL_NULL):
    p = rng.permutation(ctrl_idx); a, b = p[:N], p[N:2*N]
    d0, d1 = pair_distance(Zall, a, b)
    null_h0.append(d0); null_h1.append(d1)
null_h0 = np.array(null_h0); null_h1 = np.array(null_h1)
print(f"[control-null topo] H0 mean={null_h0.mean():.3f} p95={np.quantile(null_h0,0.95):.3f} | "
      f"H1 mean={null_h1.mean():.3f} p95={np.quantile(null_h1,0.95):.3f}  ({time.time()-t0:.0f}s)")

# ================= 3. per-gene observed + repeated sampling + permutation =================
def observed_and_repeats(g):
    # fixed observed draw
    rg = np.random.default_rng(SEED + abs(hash(g)) % 10000)
    c = rg.choice(ctrl_idx, N, replace=False)
    pcells = pert_idx[g]
    pp = rg.choice(pcells, N, replace=False) if len(pcells) >= N else pcells
    obs0, obs1 = pair_distance(Zall, c, pp)
    # repeated sampling variation
    r0, r1 = [], []
    for k in range(N_REPEAT):
        rr = np.random.default_rng(SEED + 1000*k + abs(hash(g)) % 1000)
        c2 = rr.choice(ctrl_idx, N, replace=False)
        p2 = rr.choice(pcells, N, replace=False) if len(pcells) >= N else pcells
        d0, d1 = pair_distance(Zall, c2, p2)
        r0.append(d0); r1.append(d1)
    return obs0, obs1, np.array(r0), np.array(r1), c, pp

def permutation_test(c_idx, p_idx, obs0, obs1, Z):
    pool = np.concatenate([c_idx, p_idx]); n = len(c_idx)
    perm0, perm1 = [], []
    for _ in range(N_PERM):
        q = rng.permutation(pool); a, b = q[:n], q[n:2*n]
        d0, d1 = pair_distance(Z, a, b)
        perm0.append(d0); perm1.append(d1)
    perm0 = np.array(perm0); perm1 = np.array(perm1)
    p0 = (1 + np.sum(perm0 >= obs0)) / (1 + N_PERM)
    p1 = (1 + np.sum(perm1 >= obs1)) / (1 + N_PERM)
    return float(p0), float(p1)

rows = []
for g in pilot_ids:
    obs0, obs1, r0, r1, c_used, p_used = observed_and_repeats(g)
    p0, p1 = permutation_test(c_used, p_used, obs0, obs1, Zall)
    rows.append({
        "gene_id": g, "gene_symbol": SYM[g],
        "wass_H0": obs0, "wass_H1": obs1,
        "H0_repeat_mean": float(r0.mean()), "H0_repeat_std": float(r0.std()),
        "H1_repeat_mean": float(r1.mean()), "H1_repeat_std": float(r1.std()),
        "ctrl_null_H0_mean": float(null_h0.mean()), "ctrl_null_H0_p95": float(np.quantile(null_h0, 0.95)),
        "ctrl_null_H1_mean": float(null_h1.mean()), "ctrl_null_H1_p95": float(np.quantile(null_h1, 0.95)),
        "H0_over_null": obs0 / (null_h0.mean() + 1e-12),
        "H1_over_null": obs1 / (null_h1.mean() + 1e-12),
        "perm_p_H0": p0, "perm_p_H1": p1,
    })
    print(f"  {SYM[g]:7s} H0={obs0:.3f}(p={p0:.3f}) H1={obs1:.3f}(p={p1:.3f})")

res = pd.DataFrame(rows)

# BH-FDR across the 8 genes, per dimension
def bh(pvals):
    p = np.asarray(pvals); m = len(p); order = np.argsort(p)
    adj = np.empty(m); prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = m - rank
        prev = min(prev, p[i] * m / k); adj[i] = prev
    return np.clip(adj, 0, 1)
res["perm_q_H0"] = bh(res["perm_p_H0"]); res["perm_q_H1"] = bh(res["perm_p_H1"])

# merge centroid info from pilot effect summary
eff = pd.read_csv(RESULTS / "pilot/pilot_pca_effect_summary.csv")
res = res.merge(eff[["gene_id", "centroid_shift_l2", "dispersion_ratio_to_control", "n_rpe1_available"]],
                on="gene_id", how="left")
res["centroid_exceeds_150v150_null_p95"] = res["centroid_shift_l2"] > cs_p95
res = res.sort_values("wass_H1", ascending=False)
res.to_csv(OUT / "topological_effect_size_table.csv", index=False)
pd.DataFrame([{"centroid_null_mean": cs_mean, "centroid_null_p95": cs_p95,
               "topo_null_H0_mean": float(null_h0.mean()), "topo_null_H0_p95": float(np.quantile(null_h0, .95)),
               "topo_null_H1_mean": float(null_h1.mean()), "topo_null_H1_p95": float(np.quantile(null_h1, .95)),
               "N_per_group": N, "N_pcs": N_PCS, "n_perm": N_PERM}]
             ).to_csv(OUT / "null_baselines.csv", index=False)

# ================= 4. persistence-diagram figures (control vs each gene) =================
rg = np.random.default_rng(SEED)
c_fig = rg.choice(ctrl_idx, N, replace=False)
c0, c1 = dgms_of(Zall[c_fig])
for g in pilot_ids:
    pc = pert_idx[g]
    pp = rg.choice(pc, N, replace=False) if len(pc) >= N else pc
    g0, g1 = dgms_of(Zall[pp])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    plot_diagrams([c0, c1], ax=ax[0], show=False); ax[0].set_title("Control (non-targeting)")
    plot_diagrams([g0, g1], ax=ax[1], show=False); ax[1].set_title(f"{SYM[g]} ({g})")
    fig.suptitle(f"Persistence diagrams — control vs {SYM[g]}", y=1.02)
    fig.tight_layout(); fig.savefig(DIAG_DIR / f"pd_{SYM[g]}_{g}.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# ================= 5. comparison figure =================
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for a, (xcol, xlab) in zip(ax, [("centroid_shift_l2", "centroid shift (L2, PCA)"),
                                ("dispersion_ratio_to_control", "dispersion ratio to control")]):
    a.scatter(res[xcol], res["wass_H1"], s=60, c="#2b6cb0", zorder=3)
    for _, row in res.iterrows():
        a.annotate(row["gene_symbol"], (row[xcol], row["wass_H1"]), fontsize=8,
                   xytext=(4, 4), textcoords="offset points")
    a.axhline(np.quantile(null_h1, 0.95), ls="--", c="grey", lw=1, label="H1 control null p95")
    a.set_xlabel(xlab); a.set_ylabel("Topological effect (H1 Wasserstein)"); a.legend(fontsize=8)
fig.suptitle("Topological effect size vs conventional effect measures (K562 pilot)")
fig.tight_layout(); fig.savefig(OUT / "topo_vs_conventional.png", dpi=180, bbox_inches="tight"); plt.close(fig)

# ================= 6. robustness: drop the targeted gene column, refit PCA =================
rob_rows = []
for g in pilot_ids:
    cols = [i for i in range(len(gene_cols)) if i != gidx.get(g, -1)]
    pca_r = PCA(n_components=N_PCS, svd_solver="randomized", random_state=SEED)
    Zr = pca_r.fit_transform(Xall[:, cols]).astype(np.float32)
    rg2 = np.random.default_rng(SEED + abs(hash(g)) % 10000)
    c = rg2.choice(ctrl_idx, N, replace=False)
    pc = pert_idx[g]; pp = rg2.choice(pc, N, replace=False) if len(pc) >= N else pc
    o0, o1 = pair_distance(Zr, c, pp)
    p0, p1 = permutation_test(c, pp, o0, o1, Zr)
    rob_rows.append({"gene_id": g, "gene_symbol": SYM[g],
                     "wass_H0_drop_target": o0, "wass_H1_drop_target": o1,
                     "perm_p_H0_drop": p0, "perm_p_H1_drop": p1})
rob = pd.DataFrame(rob_rows)
rob = rob.merge(res[["gene_id", "wass_H0", "wass_H1"]], on="gene_id", how="left")
rob["H1_retained_frac"] = rob["wass_H1_drop_target"] / (rob["wass_H1"] + 1e-12)
rob.to_csv(OUT / "robustness_drop_target_gene.csv", index=False)

# ================= matched pairs summary =================
pairs = [("ENSG00000168827", "ENSG00000099821"),   # GFM1 vs POLRMT
         ("ENSG00000122026", "ENSG00000114503"),   # RPL21 vs NCBP2
         ("ENSG00000119203", "ENSG00000146007")]   # CPSF3 vs ZMAT2
ri = res.set_index("gene_id")
pair_rows = []
for a, b in pairs:
    for g in (a, b):
        row = ri.loc[g]
        pair_rows.append({"pair": f"{SYM[a]} vs {SYM[b]}", "gene": SYM[g],
                          "centroid_shift": round(row["centroid_shift_l2"], 3),
                          "wass_H0": round(row["wass_H0"], 3), "wass_H1": round(row["wass_H1"], 3),
                          "H1_over_null": round(row["H1_over_null"], 2),
                          "perm_p_H1": row["perm_p_H1"], "perm_q_H1": round(row["perm_q_H1"], 3)})
pairdf = pd.DataFrame(pair_rows)
pairdf.to_csv(OUT / "matched_pair_comparison.csv", index=False)

print("\n===== TOPOLOGICAL EFFECT TABLE =====")
print(res[["gene_symbol", "centroid_shift_l2", "wass_H0", "wass_H1", "H1_over_null",
           "perm_p_H1", "perm_q_H1"]].round(3).to_string(index=False))
print("\n===== MATCHED PAIRS =====")
print(pairdf.to_string(index=False))
print("\n===== ROBUSTNESS (drop targeted gene) =====")
print(rob[["gene_symbol", "wass_H1", "wass_H1_drop_target", "H1_retained_frac", "perm_p_H1_drop"]].round(3).to_string(index=False))
print("\nOutputs in", OUT)
for p in sorted(OUT.rglob("*")):
    if p.is_file():
        print("  ", p.relative_to(OUT))
