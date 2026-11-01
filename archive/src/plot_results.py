"""Heatmaps: (a) synthetic detection rates with known truth, (b) real knockdowns."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from methods import METHODS, LABELS
from synthetic import ORDER, SIGMAS, D

RESULTS = Path(__file__).resolve().parents[1] / "results"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
COND_LABEL = {"null": "Null (same Gaussian)", "translation": "Translation",
              "scaling": "Scaling ×1.3", "anisotropy": "Anisotropy (PR 16.6→10)",
              "heavy_tails": "Heavy tails (t, df 5)",
              "partial_kd": "Partial knockdown (70% shifted)",
              "ring_k562": "Loop in top-2 PCs (26% var.)", "bimodal_k562": "Two clusters on PC1",
              **{f"ring_s{s}": f"Loop + noise σ={s} ({(1 + 2*s*s)/(1 + D*s*s):.0%} var.)" for s in SIGMAS}}


def heat(ax, M, rows, title, marks=None, sep=None):
    ax.imshow(M, cmap="Greys", vmin=0, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            t = f"{v * 100:.0f}" + (" *" if marks is not None and marks[i, j] else "")
            ax.text(j, i, t, ha="center", va="center", fontsize=7, color="white" if v > 0.55 else "black")
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels([LABELS[m] for m in METHODS], rotation=40, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows)
    ax.set_xticks(np.arange(-.5, len(METHODS)), minor=True)
    ax.set_yticks(np.arange(-.5, len(rows)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    for x in (1.5, 6.5):
        ax.axvline(x, color="black", linewidth=1)
    for y in sep or []:
        ax.axhline(y, color="black", linewidth=1)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold", pad=24)


def main(results):
    syn = pd.read_csv(results / "synthetic.csv", keep_default_na=False)
    conds = ORDER
    Ms = syn.pivot(index="condition", columns="method", values="detection_rate").loc[conds, METHODS].values
    real = pd.read_csv(results / "real.csv")
    order = list(dict.fromkeys(zip(real.cell_line, real.group)))
    Mr = real.pivot_table(index=["cell_line", "group"], columns="method", values="detect_frac", sort=False).loc[order, METHODS].values
    Qr = real.pivot_table(index=["cell_line", "group"], columns="method", values="q_primary", sort=False).reindex(order)[METHODS].values
    rlab = [f"{cl}  {g}" for cl, g in order]
    nk = sum(cl == "K562" for cl, _ in order)

    fig, (a, b) = plt.subplots(2, 1, figsize=(8.2, 9.6), gridspec_kw={"height_ratios": [len(conds), len(order)]})
    heat(a, Ms, [COND_LABEL[c] for c in conds],
         "a  Synthetic clouds with known truth: % of 40 replicates with p < 0.05", sep=[4.5, 5.5])
    for y, t in ((2, "geometry only\n(should stay ~5)"), (5, "mixture"), (8.5, "topology changed\n(should be high)")):
        a.text(1.02, y, t, transform=a.get_yaxis_transform(), ha="left", va="center", fontsize=7, style="italic")
    heat(b, Mr, rlab, "b  Real knockdowns: % of 15 draws with p < 0.05   (* primary draw BH q < 0.05)",
         marks=np.nan_to_num(Qr, nan=1) < 0.05, sep=[nk - 0.5])
    for ax in (a, b):
        for x, t in ((0.5, "conventional"), (4, "persistent homology"), (8, "covariance-matched\nsurrogate")):
            ax.text(x, -0.6, t, ha="center", va="bottom", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(results / "benchmark.png", dpi=220, bbox_inches="tight")
    print("saved", results / "benchmark.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    args = parser.parse_args()
    main(args.results)
