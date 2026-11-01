"""Step 2b -- the same battery on the real K562 and RPE1 knockdowns.

For each perturbation we draw ``R`` times: N control cells vs N perturbed cells,
and run every test. Draw 0 is the primary analysis (BH-corrected across genes
within a cell line); the fraction of draws that reject shows how robust a hit is.
A control-vs-control group (controls split in two) is the on-data sanity check.

Beyond the p-values we save the interpretable loop readout (Step 3): the leading
H1 lifetime for control vs perturbed, its signed change, and -- for the primary
draw -- a bootstrap interval and the persistence diagrams themselves.

    python experiments/03_real.py            # 15 draws per group
    python experiments/03_real.py 5          # quicker
"""
import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import config
from topoperturb import data, features, homology, plotting, twosample

CONTROL = data.CONTROL
N, SEED = config.N_CELLS, config.SEED
CLOUDS = {cl: data.load_clouds(config.cloud_path(cl))[0] for cl in ("k562", "rpe1")}
# (cell_line, group) in display order; a control-vs-control null per line
GROUPS = [(cl, g) for cl in CLOUDS for g in config.TARGETS[cl] + ["control_vs_control"]]


def sample(cell_line, group, draw):
    """Draw a control and a perturbed cloud for one (group, draw); reproducible.

    Controls come only from the testing pool: the reference controls that defined
    the PCA axes are excluded, so no cell both built and is judged in the space.
    """
    Z, labels, _rows, reference = CLOUDS[cell_line]
    control_idx = data.testing_controls(labels, reference)
    rng = np.random.default_rng([SEED, 77, GROUPS.index((cell_line, group)), draw])
    if group == "control_vs_control":
        pick = rng.choice(control_idx, 2 * N, replace=False)
        return Z[pick[:N]], Z[pick[N:]], rng
    control = rng.choice(control_idx, N, replace=False)
    perturbed = rng.choice(np.where(labels == group)[0], N, replace=False)
    return Z[control], Z[perturbed], rng


def run(task):
    """One draw of one group: battery p-values plus the signed loop change."""
    cell_line, group, draw = task
    control, perturbed, rng = sample(cell_line, group, draw)
    row = {"cell_line": cell_line, "group": config.symbol(group), "draw": draw,
           **twosample.compare(control, perturbed, rng)}
    row.update(features.signed_change(control, perturbed))
    return row


def benjamini_hochberg(pvalues):
    """BH-adjusted q-values for a set of p-values."""
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def summarise(draws):
    """Per group and method: detection fraction, median p, primary p, BH q."""
    rows = []
    for (cell_line, group), subset in draws.groupby(["cell_line", "group"], sort=False):
        for method in twosample.METHODS:
            p = subset[f"p_{method}"]
            rows.append({"cell_line": cell_line, "group": group, "method": method,
                         "detect_frac": (p < 0.05).mean(), "median_p": p.median(),
                         "p_primary": subset.loc[subset.draw == 0, f"p_{method}"].iloc[0]})
    table = pd.DataFrame(rows)
    is_gene = ~table.group.str.startswith(("control", "pos."))
    for (_cl, _m), idx in table[is_gene].groupby(["cell_line", "method"]).groups.items():
        table.loc[idx, "q_primary"] = benjamini_hochberg(table.loc[idx, "p_primary"])
    return table


def feature_report(out):
    """Primary-draw loop readout + saved diagrams for the real perturbations."""
    (out / "diagrams").mkdir(parents=True, exist_ok=True)
    rows = []
    for cell_line, group in GROUPS:
        if group == "control_vs_control":
            continue
        control, perturbed, rng = sample(cell_line, group, 0)
        row = {"cell_line": cell_line, "group": config.symbol(group)}
        row.update(features.signed_change(control, perturbed))
        row.update(features.bootstrap_leading_h1(perturbed, rng=rng))
        rows.append(row)
        h0c, h1c = homology.summarize(control, ("diffusion",))["diffusion"].values()
        h0p, h1p = homology.summarize(perturbed, ("diffusion",))["diffusion"].values()
        fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=(7, 3.4), sharex=True, sharey=True)
        plotting.persistence_diagram(ax1, h0c, h1c, f"{cell_line} control")
        plotting.persistence_diagram(ax2, h0p, h1p, f"{cell_line} {config.symbol(group)}")
        plotting.save(fig, out / "diagrams" / f"{cell_line}_{config.symbol(group)}.png")
    pd.DataFrame(rows).to_csv(out / "real_features.csv", index=False)


def main(repeats, out):
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(cl, g, d) for d in range(repeats) for cl, g in GROUPS]
    with Pool(8) as pool:
        draws = pd.DataFrame(pool.map(run, tasks, chunksize=1))
    draws.to_csv(out / "real_draws.csv", index=False)
    table = summarise(draws)
    table.to_csv(out / "real.csv", index=False)
    feature_report(out)

    order = list(dict.fromkeys(zip(table.cell_line, table.group)))
    detect = table.pivot_table(index=["cell_line", "group"], columns="method",
                               values="detect_frac", sort=False).reindex(order)[twosample.METHODS]
    q = table.pivot_table(index=["cell_line", "group"], columns="method",
                          values="q_primary", sort=False).reindex(order)[twosample.METHODS]
    fig, ax = plotting.plt.subplots(figsize=(7.5, 5.5))
    plotting.heatmap(ax, detect.values, [f"{cl}  {g}" for cl, g in order], twosample.METHODS,
                     "Real knockdowns: % of draws with p < 0.05   (* BH q < 0.05 on primary draw)",
                     marks=np.nan_to_num(q.values, nan=1) < 0.05)
    plotting.save(fig, out / "real.png")
    print(detect.round(2).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repeats", type=int, nargs="?", default=15)
    parser.add_argument("--output", type=Path, default=config.RESULTS)
    args = parser.parse_args()
    t = time.time()
    main(args.repeats, args.output)
    print(f"{time.time() - t:.0f}s")
