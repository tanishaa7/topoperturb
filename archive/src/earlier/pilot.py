"""Reproduce the count-selected pilot, balanced PCA and optional Mapper views.

Uses the recovered expression caches plus the saved full count inventory.
The historical centroid null splits 150 controls into 75-vs-75 halves; its
n_per_group output labels the pilot size, not the null comparison size.
"""

if __name__ == "__main__":
    from common import configure_run
    args = configure_run(__doc__, mapper=True)

from pathlib import Path
from common import DATA, RESULTS, load, load_cache
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

SEED = 42
OUT = RESULTS / "pilot"
OUT.mkdir(exist_ok=True)

# Recovered caches preserve expression rows and measured gene names. The full
# inventories include excluded cells that are deliberately absent from caches.
Xk, lk, gk, _ = load()
_, _, gr = load_cache(DATA / "expression_rpe1.npz")
gene_cols = {"K562": gk, "RPE1": gr}
inventory = pd.read_csv(DATA / "cell_counts.csv")
counts = {}
for line in ("K562", "RPE1"):
    table = inventory[inventory.cell_line == line].drop(columns="cell_line")
    counts[line] = table.set_index("intervention").n_cells
    table.to_csv(OUT / f"{line.lower()}_perturbation_counts.csv", index=False)
    print(f"{line}: {int(counts[line].sum()):,} inventoried cells; {len(gene_cols[line])} genes")

# ---- 3. shared targets + balanced pilot ----
def target_table(line):
    genes = set(gene_cols[line])
    s = counts[line][counts[line].index.isin(genes)]
    return pd.DataFrame({"gene_id": s.index.astype(str), f"n_{line.lower()}": s.to_numpy(dtype=int)})

k = target_table("K562"); r = target_table("RPE1")
shared = k.merge(r, on="gene_id", how="inner").sort_values("n_k562", ascending=False)
shared.to_csv(OUT / "shared_perturbations.csv", index=False)
print(f"\nShared measured-and-targeted genes: {len(shared)}")

def spread_pick(df, n, count_col="n_k562"):
    df = df.sort_values([count_col, "gene_id"]).reset_index(drop=True)
    if df.empty or n <= 0:
        return df.iloc[:0].copy()
    n = min(n, len(df))
    ranks = np.linspace(0, len(df) - 1, n).round().astype(int)
    return df.iloc[np.unique(ranks)].copy()

shared_candidates = shared[(shared.n_k562 >= 200) & (shared.n_rpe1 >= 150)]
if len(shared_candidates) < 4:
    shared_candidates = shared[(shared.n_k562 >= 150) & (shared.n_rpe1 >= 100)]
k_only = k[~k.gene_id.isin(set(shared.gene_id))]
k_only_candidates = k_only[k_only.n_k562 >= 200]
if len(k_only_candidates) < 4:
    k_only_candidates = k_only[k_only.n_k562 >= 150]

pilot_shared = spread_pick(shared_candidates, min(4, len(shared_candidates)))
pilot_only = spread_pick(k_only_candidates, min(4, len(k_only_candidates)))
pilot_shared["group"] = "shared_RPE1"
pilot_only["group"] = "K562_only"
pilot = pd.concat([pilot_shared, pilot_only], ignore_index=True)
if len(pilot) < 8:
    remaining = k[~k.gene_id.isin(set(pilot.gene_id)) & (k.n_k562 >= 100)]
    filler = spread_pick(remaining, 8 - len(pilot))
    filler["group"] = "additional_K562"
    pilot = pd.concat([pilot, filler], ignore_index=True)
pilot = pilot.merge(r, on="gene_id", how="left", suffixes=("", "_lookup"))
if "n_rpe1_lookup" in pilot.columns:
    pilot["n_rpe1"] = pilot.get("n_rpe1", pilot["n_rpe1_lookup"]).fillna(pilot["n_rpe1_lookup"])
    pilot = pilot.drop(columns=["n_rpe1_lookup"])
pilot["n_rpe1"] = pilot["n_rpe1"].fillna(0).astype(int)
pilot = pilot[["gene_id", "n_k562", "n_rpe1", "group"]]
assert pilot.gene_id.is_unique
# Reuse the supplied symbols; do not invent an annotation for new targets.
annotation = pd.read_csv(DATA / "pilot.csv").set_index("gene_id").gene_symbol
pilot["gene_symbol"] = pilot.gene_id.map(annotation)
if pilot.gene_symbol.isna().any():
    raise ValueError("Pilot selection changed; supply a versioned symbol annotation before proceeding")
pilot.to_csv(OUT / "pilot_perturbations.csv", index=False)
print("\nPilot shortlist:")
print(pilot.to_string(index=False))

# ---- 4. load control + pilot cells ----
pilot_ids = pilot.gene_id.tolist()
if not set(pilot_ids).issubset(set(lk)):
    raise ValueError("Selected pilot is not covered by the recovered K562 cache")
mask = np.isin(lk, ["non-targeting", *pilot_ids])
pilot_data = pd.DataFrame(Xk[mask], columns=gene_cols["K562"])
pilot_data["intervention"] = lk[mask]

# ---- 5. shared PCA ----
N_PER_GROUP = min(150, int(pilot_data.intervention.value_counts().loc[pilot_ids].min()))
assert N_PER_GROUP >= 30
sampled = []
for j, label in enumerate(["non-targeting", *pilot_ids]):
    sub = pilot_data.loc[pilot_data.intervention == label]
    sampled.append(sub.sample(n=N_PER_GROUP, random_state=SEED + j, replace=False))
balanced = pd.concat(sampled, ignore_index=True)
X = balanced.loc[:, gene_cols["K562"]].to_numpy(dtype=np.float32, copy=True)
assert np.isfinite(X).all()
Y = balanced.intervention.astype(str).to_numpy()
n_components = min(30, X.shape[1], X.shape[0] - 1)
pca = PCA(n_components=n_components, svd_solver="randomized", random_state=SEED)
Z = pca.fit_transform(X)

plt.figure(figsize=(8.5, 6))
focus = ["non-targeting", *pilot_ids[:3]]
for label in focus:
    loc = Y == label
    plt.scatter(Z[loc, 0], Z[loc, 1], s=12, alpha=0.65, label=label)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
plt.title(f"CausalBench K562 pilot: balanced PCA, n={N_PER_GROUP}/group")
plt.legend(title="Intervention", fontsize=8)
plt.tight_layout()
plt.savefig(OUT / "pilot_pca.png", dpi=220)
print("\nPCA shape:", Z.shape, "ExplVar first10:", round(float(pca.explained_variance_ratio_[:10].sum()), 3))

# ---- 6. effect summary WITH control-vs-control null baseline (ADDED) ----
z_ctrl = Z[Y == "non-targeting"]
c_ctrl = z_ctrl.mean(axis=0)
s_ctrl = np.mean(np.linalg.norm(z_ctrl - c_ctrl, axis=1))

# ADDED: control-vs-control baseline. Repeatedly split controls into two
# disjoint halves and measure the centroid shift between halves -> the
# centroid displacement expected from sampling noise alone.
rng = np.random.default_rng(SEED)
n_half = len(z_ctrl) // 2
null_shifts = []
for _ in range(500):
    perm = rng.permutation(len(z_ctrl))
    a, b = z_ctrl[perm[:n_half]], z_ctrl[perm[n_half:2 * n_half]]
    null_shifts.append(np.linalg.norm(a.mean(0) - b.mean(0)))
null_shifts = np.array(null_shifts)
null_mean = float(null_shifts.mean())
null_p95 = float(np.quantile(null_shifts, 0.95))
print(f"\nControl-vs-control null centroid shift: mean={null_mean:.3f}, 95th pct={null_p95:.3f}")

rows = []
for gene in pilot_ids:
    z = Z[Y == gene]
    c = z.mean(axis=0)
    shift = float(np.linalg.norm(c - c_ctrl))
    rows.append({"gene_id": gene,
                 "n_pilot_cells": len(z),
                 "centroid_shift_l2": shift,
                 "shift_over_null_mean": shift / (null_mean + 1e-12),
                 "exceeds_null_p95": shift > null_p95,
                 "mean_radial_dispersion": float(np.mean(np.linalg.norm(z - c, axis=1))),
                 "dispersion_ratio_to_control": float(np.mean(np.linalg.norm(z - c, axis=1)) / (s_ctrl + 1e-12)),
                 "n_rpe1_available": int(pilot.set_index("gene_id").loc[gene, "n_rpe1"])})
effects = pd.DataFrame(rows).sort_values("centroid_shift_l2", ascending=False)
effects.attrs["null_mean"] = null_mean
effects.to_csv(OUT / "pilot_pca_effect_summary.csv", index=False)
# store the null row alongside as a small companion file
pd.DataFrame([{"null_centroid_shift_mean": null_mean, "null_centroid_shift_p95": null_p95,
               "control_mean_radial_dispersion": float(s_ctrl), "n_per_group": N_PER_GROUP,
               "null_cells_per_group": n_half}]
             ).to_csv(OUT / "pilot_control_null_baseline.csv", index=False)
print("\nEffect summary (centroid shift vs control, with null reference):")
print(effects.round(3).to_string(index=False))

# ---- 7. optional Mapper ----
if not args.skip_mapper:
    import kmapper as km
    mapper = km.KeplerMapper(verbose=0)
    HTML_GENES = pilot_ids[:3]
    for gene in HTML_GENES:
        mask = np.isin(Y, ["non-targeting", gene])
        z = Z[mask]; y = Y[mask]
        embed = z[:, :min(10, z.shape[1])]
        nn = NearestNeighbors(n_neighbors=min(8, len(embed) - 1)).fit(embed)
        dists, _ = nn.kneighbors(embed)
        eps = max(float(np.quantile(dists[:, -1], 0.80)) * 1.25, 1e-6)
        graph = mapper.map(lens=z[:, :2], X=embed,
                           cover=km.Cover(n_cubes=8, perc_overlap=0.4),
                           clusterer=DBSCAN(eps=eps, min_samples=4))
        output_path = OUT / f"mapper_{gene}.html"
        mapper.visualize(graph, path_html=str(output_path),
                         title=f"K562 control versus {gene}",
                         color_values=(y == gene).astype(float),
                         color_function_name=f"Perturbed ({gene}) fraction")
        print(f"{gene}: {len(graph['nodes'])} Mapper nodes; eps={eps:.3f}; saved {output_path.name}")

print("\nAll outputs in:", OUT)
for p in sorted(OUT.iterdir()):
    print("  ", p.name)
