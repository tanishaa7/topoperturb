"""Step 2c -- does each test keep its promised false-positive rate?

A test at p < 0.05 should reject about 5% of the time when nothing differs. We
check this on real cells by comparing controls with controls: for each cell line
and several group sizes we draw many control-vs-control pairs and record how often
each method rejects, with a Wilson interval, and whether its p-values look uniform
(a KS test against Uniform(0,1)). Methods whose rate sits well above 5% are not
trustworthy as evidence, however strong they look on the real knockdowns.

Cells come from the testing controls only: the reference controls that fitted the
PCA are held out, so the representation is independent of the comparison and this
is a true null rather than a test of cells against their own axes.

The draws reuse cells from the same control pool, so they are not independent; the
Wilson interval reported here therefore *understates* the true uncertainty and
should be read as a lower bound on it.

    python experiments/04_calibration.py              # 200 draws, sizes 75/100/150
    python experiments/04_calibration.py --draws 50
    python experiments/04_calibration.py --n-pcs 50   # representation sensitivity
"""
import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kstest

import config
from topoperturb import data, plotting, twosample

CONTROL = data.CONTROL
SEED = config.SEED
SIZES = [75, 100, 150]
CELL_LINES = ("k562", "rpe1")
_CLOUDS = {}


def clouds(cell_line, n_pcs):
    """Clouds for one cell line and PC count, cached per worker process."""
    if (cell_line, n_pcs) not in _CLOUDS:
        _CLOUDS[cell_line, n_pcs] = data.load_clouds(config.cloud_path(cell_line, n_pcs))[0]
    return _CLOUDS[cell_line, n_pcs]


def run(task):
    """One control-vs-control draw at a given size: all methods' p-values.

    Cells are drawn from the testing controls only -- the PCA reference controls are
    excluded, so this really is a null: the axes were fitted on cells that cannot
    appear on either side of the comparison. The draw does not depend on ``n_pcs``,
    so different representations are compared on identical cells and permutations.
    """
    cell_line, size, draw, n_pcs = task
    Z, labels, _rows, reference = clouds(cell_line, n_pcs)
    control_idx = data.testing_controls(labels, reference)
    rng = np.random.default_rng([SEED, 13, SIZES.index(size), draw])
    pick = rng.choice(control_idx, 2 * size, replace=False)
    row = {"cell_line": cell_line, "size": size, "draw": draw, "n_pcs": n_pcs}
    for method, (p, _obs) in twosample.permutation_test(
            twosample.all_statistics, Z[pick[:size]], Z[pick[size:]], rng=rng).items():
        row[f"p_{method}"] = p
    return row


def summarise(draws):
    """False-positive rate, Wilson interval and p-value uniformity per cell/size/method."""
    rows = []
    for (cell_line, size), subset in draws.groupby(["cell_line", "size"]):
        for method in twosample.METHODS:
            p = subset[f"p_{method}"].to_numpy()
            k = int((p < 0.05).sum())
            lo, hi = twosample.wilson_interval(k, len(p))
            rows.append({"cell_line": cell_line, "size": size, "method": method,
                         "n_pcs": int(subset.n_pcs.iloc[0]),
                         "false_positive_rate": k / len(p), "ci_low": lo, "ci_high": hi,
                         "n_draws": len(p), "uniform_ks_p": float(kstest(p, "uniform").pvalue)})
    return pd.DataFrame(rows)


def main(draws_per_cell, out, n_pcs):
    out.mkdir(parents=True, exist_ok=True)
    tag = "" if n_pcs == config.N_PCS else f"_{n_pcs}pc"
    tasks = [(cl, size, d, n_pcs)
             for cl in CELL_LINES for size in SIZES for d in range(draws_per_cell)]
    with Pool(8) as pool:
        draws = pd.DataFrame(pool.map(run, tasks, chunksize=1))
    draws.to_csv(out / f"calibration_draws{tag}.csv", index=False)
    table = summarise(draws)
    table.to_csv(out / f"calibration{tag}.csv", index=False)

    fig, axes = plotting.plt.subplots(1, len(CELL_LINES), figsize=(5 * len(CELL_LINES), 3.6),
                                      squeeze=False, sharey=True)
    for ax, cell_line in zip(axes[0], CELL_LINES):
        for method in twosample.METHODS:
            sub = table[(table.cell_line == cell_line) & (table.method == method)].sort_values("size")
            plotting.calibration_curve(ax, sub["size"], sub.false_positive_rate,
                                       sub.ci_low, sub.ci_high, method)
        ax.axhline(0.05, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(f"{cell_line}: control-vs-control false-positive rate "
                     f"({n_pcs} PCs)", fontsize=9)
    axes[0][-1].legend(fontsize=6, loc="upper right")
    plotting.save(fig, out / f"calibration{tag}.png")
    for value, caption in [("false_positive_rate", "false-positive rate at p < 0.05"),
                           ("uniform_ks_p", "KS p-value vs Uniform(0,1)")]:
        print(f"\n{caption} ({n_pcs} PCs, {draws_per_cell} draws per size):")
        print(table.pivot_table(index=["cell_line", "size"], columns="method",
                                values=value).round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--draws", type=int, default=200, help="control-vs-control pairs per size")
    parser.add_argument("--output", type=Path, default=config.RESULTS)
    parser.add_argument("--n-pcs", type=int, default=config.N_PCS,
                        help=f"which prepared representation to test (default {config.N_PCS})")
    args = parser.parse_args()
    t = time.time()
    main(args.draws, args.output, args.n_pcs)
    print(f"{time.time() - t:.0f}s")
