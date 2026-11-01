"""Small, consistent figure helpers. Each experiment saves one figure with these."""
import matplotlib

matplotlib.use("Agg")       # file output only; no interactive display needed
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})


def heatmap(ax, matrix, row_labels, col_labels, title, marks=None):
    """Grey detection-rate heatmap in [0, 1]; ``marks`` adds a '*' per cell."""
    ax.imshow(matrix, cmap="Greys", vmin=0, vmax=1, aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text = f"{value * 100:.0f}" + (" *" if marks is not None and marks[i, j] else "")
            ax.text(j, i, text, ha="center", va="center", fontsize=7,
                    color="white" if value > 0.55 else "black")
    ax.set_xticks(range(len(col_labels)), col_labels, rotation=40, ha="right")
    ax.set_yticks(range(len(row_labels)), row_labels)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold")


def persistence_diagram(ax, h0, h1, title=""):
    """Scatter a persistence diagram with the birth = death line for reference."""
    high = 1.0
    for dgm in (h0, h1):
        if len(dgm):
            high = max(high, float(dgm[:, 1].max()))
    ax.plot([0, high], [0, high], color="0.7", linewidth=0.8)
    if len(h0):
        ax.scatter(h0[:, 0], h0[:, 1], s=12, color="#4477aa", label="H0")
    if len(h1):
        ax.scatter(h1[:, 0], h1[:, 1], s=16, color="#cc6677", marker="^", label="H1")
    ax.set_xlabel("birth")
    ax.set_ylabel("death")
    ax.set_title(title, fontsize=8)
    ax.legend(fontsize=7, loc="lower right")


def calibration_curve(ax, sizes, rates, lows, highs, label):
    """Rejection rate vs sample size with a shaded confidence band for one method."""
    sizes, rates = np.asarray(sizes), np.asarray(rates)
    line, = ax.plot(sizes, rates, marker="o", label=label)
    ax.fill_between(sizes, lows, highs, color=line.get_color(), alpha=0.15)
    ax.set_xlabel("cells per group")
    ax.set_ylabel("false-positive rate")


def save(fig, path):
    """Tight-layout and write a figure, printing where it went."""
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    print("saved", path)
