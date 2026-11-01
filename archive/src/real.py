"""The same test battery applied to the real K562 and RPE1 knockdowns.

Per group, R draws of 150 control vs 150 perturbed cells (30 PCs). Draw 0 is
the primary analysis (BH across genes within each cell line); the fraction of
draws with p < 0.05 shows robustness. Draws reuse cells, so they are not
independent replicates. Control-vs-control groups (300 controls split in two)
give the empirical false-positive rate on real cells.
"""
import argparse
import json
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
from methods import battery, METHODS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
META = json.loads((ROOT / "data/metadata.json").read_text())
_z = np.load(ROOT / "data/clouds.npz")
CLOUDS = {"K562": (_z["Zk"], _z["lk"].astype(str)), "RPE1": (_z["Zr"], _z["lr"].astype(str))}
SYMS = {"K562": META["k562_symbols"], "RPE1": META["rpe1_symbols"]}
N, SEED = 150, 42
GROUPS = [(cl, g) for cl in CLOUDS for g in list(SYMS[cl]) + ["control_vs_control"]]


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.clip(q, 0, 1)
    return out


def run(task):
    gi, k = task
    cl, g = GROUPS[gi]
    Z, lab = CLOUDS[cl]
    ctrl = np.where(lab == "non-targeting")[0]
    rng = np.random.default_rng([SEED, 77, gi, k])
    if g == "control_vs_control":
        q = rng.choice(ctrl, 2 * N, replace=False)
        c, p = q[:N], q[N:]
    else:
        c = rng.choice(ctrl, N, replace=False)
        p = rng.choice(np.where(lab == g)[0], N, replace=False)
    name = SYMS[cl].get(g, "control vs control")
    return {"cell_line": cl, "group": name, "draw": k, **battery(Z[c], Z[p], rng)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repeats", type=int, nargs="?", default=15)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    OUT = args.output
    OUT.mkdir(parents=True, exist_ok=True)
    t = time.time()
    tasks = [(gi, k) for k in range(args.repeats) for gi in range(len(GROUPS))]
    with Pool(8) as pool:
        rows = pool.map(run, tasks, chunksize=1)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "real_draws.csv", index=False)
    out = []
    for (cl, g), s in df.groupby(["cell_line", "group"], sort=False):
        for m in METHODS:
            out.append({"cell_line": cl, "group": g, "method": m,
                        "detect_frac": (s[m] < 0.05).mean(), "median_p": s[m].median(),
                        "p_primary": s.loc[s.draw == 0, m].iloc[0]})
    tab = pd.DataFrame(out)
    genes = ~tab.group.str.startswith(("control", "pos."))
    for (cl, m), idx in tab[genes].groupby(["cell_line", "method"]).groups.items():
        tab.loc[idx, "q_primary"] = bh(tab.loc[idx, "p_primary"])
    tab.to_csv(OUT / "real.csv", index=False)
    pd.set_option("display.width", 250)
    for v in ("detect_frac", "q_primary"):
        print(f"== {v} ==")
        print(tab.pivot_table(index=["cell_line", "group"], columns="method", values=v, sort=False)
              [METHODS].round(2).to_string())
    print(f"{time.time() - t:.0f}s")
