"""Compare synthetic geometry, ring and mixture alternatives (150 cells, 30 PCs).

The base Gaussian uses the saved K562 control covariance spectrum. Conditions
change translation, scale, anisotropy, tails, rings or population mixtures.
Covariances are matched in expectation where specified, not sample by sample.
Noisy rings/mixtures describe latent or density structure; see the research audit
before interpreting the historical `true_topology_change` labels as literal truth.

Usage: python src/synthetic.py [repeats] [condition ...] --output results/local
Named conditions are rerun and merged into an existing complete results file.
"""
import argparse
import json
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from methods import battery, METHODS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
META = json.loads((ROOT / "data/metadata.json").read_text())
LAM = np.array(META["lam"])
N, D, SEED = 150, 30, 42
SIGMAS = [0.1, 0.2, 0.3, 0.5]


def pr_of(l):
    return l.sum() ** 2 / (l ** 2).sum()


def sharpen(l, target):
    g = brentq(lambda g: pr_of(l ** g) - target, 1.0, 5.0)
    out = l ** g
    return out * l.sum() / out.sum()


LAM_ANISO = sharpen(LAM, 10.0)
CONDITIONS = {
    "null": ("geometry", "none"),
    "translation": ("geometry", "none"),
    "scaling": ("geometry", "none"),
    "anisotropy": ("geometry", "none"),
    "heavy_tails": ("geometry", "none (non-Gaussian)"),
    "ring_k562": ("topology", "loop (H1)"),
    "bimodal_k562": ("topology", "two clusters (H0)"),
    **{f"ring_s{s}": ("topology", "loop (H1)") for s in SIGMAS},
    "partial_kd": ("heterogeneity", "none (mixture)"),     # appended last so seeds of the others are unchanged
}
ORDER = ["null", "translation", "scaling", "anisotropy", "heavy_tails", "partial_kd",
         "ring_k562", "bimodal_k562"] + [f"ring_s{s}" for s in SIGMAS]


def gauss(rng, lam, n=N):
    return rng.normal(size=(n, len(lam))) * np.sqrt(lam)


def make(cond, rng):
    if cond.startswith("ring_s"):
        s = float(cond[6:])
        th = rng.uniform(0, 2 * np.pi, N)
        B = rng.normal(scale=s, size=(N, D))
        B[:, 0] += np.cos(th)
        B[:, 1] += np.sin(th)
        lam = np.full(D, s ** 2)
        lam[:2] += 0.5
        return gauss(rng, lam), B
    A = gauss(rng, LAM)
    if cond == "null":
        B = gauss(rng, LAM)
    elif cond == "translation":
        u = rng.normal(size=D)
        B = gauss(rng, LAM) + META["median_centroid_shift"] * u / np.linalg.norm(u)
    elif cond == "scaling":
        B = 1.3 * gauss(rng, LAM)
    elif cond == "anisotropy":
        B = gauss(rng, LAM_ANISO)
    elif cond == "heavy_tails":
        df = 5
        w = rng.chisquare(df, N) / df
        B = gauss(rng, LAM) * np.sqrt((df - 2) / df) / np.sqrt(w)[:, None]
    elif cond == "ring_k562":
        th = rng.uniform(0, 2 * np.pi, N)
        B = gauss(rng, LAM)
        B[:, 0] = np.sqrt(2 * LAM[0]) * np.cos(th)
        B[:, 1] = np.sqrt(2 * LAM[1]) * np.sin(th)
    elif cond == "bimodal_k562":
        d = 0.9 * np.sqrt(LAM[0])
        B = gauss(rng, LAM)
        B[:, 0] = rng.choice([-d, d], N) + rng.normal(scale=np.sqrt(LAM[0] - d ** 2), size=N)
    elif cond == "partial_kd":
        u = rng.normal(size=D)
        B = gauss(rng, LAM)
        B[: int(0.7 * N)] += META["median_centroid_shift"] * u / np.linalg.norm(u)
    return A, B


def run(task):
    ci, cond, rep = task
    rng = np.random.default_rng([SEED, ci, rep])
    A, B = make(cond, rng)
    return {"condition": cond, "rep": rep, **battery(A, B, rng)}


def wilson(k, n, z=1.96):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def summarise(df):
    rows = []
    for cond in ORDER:
        fam, truth = CONDITIONS[cond]
        s = df[df.condition == cond]
        for m in METHODS:
            k = int((s[m] < 0.05).sum())
            lo, hi = wilson(k, len(s))
            rows.append({"condition": cond, "family": fam, "true_topology_change": truth,
                         "method": m, "detection_rate": k / len(s), "ci_low": lo, "ci_high": hi,
                         "n_reps": len(s), "median_p": s[m].median()})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repeats", type=int, nargs="?", default=40)
    parser.add_argument("conditions", nargs="*")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    if set(args.conditions) - set(CONDITIONS):
        parser.error("unknown condition; choose from " + ", ".join(CONDITIONS))
    OUT = args.output
    OUT.mkdir(parents=True, exist_ok=True)
    if args.conditions and not (OUT / "synthetic_draws.csv").exists():
        parser.error("a named rerun needs an existing complete synthetic_draws.csv in --output")
    t = time.time()
    only = args.conditions or list(CONDITIONS)
    tasks = [(i, c, r) for i, c in enumerate(CONDITIONS) if c in only for r in range(args.repeats)]
    with Pool(8) as pool:
        rows = pool.map(run, tasks, chunksize=1)
    df = pd.DataFrame(rows)
    if args.conditions:
        old = pd.read_csv(OUT / "synthetic_draws.csv", keep_default_na=False)   # keep "null" as text
        df = pd.concat([old[~old.condition.isin(only)], df], ignore_index=True)
    df.to_csv(OUT / "synthetic_draws.csv", index=False)
    tab = summarise(df)
    tab.to_csv(OUT / "synthetic.csv", index=False)
    pd.set_option("display.width", 250)
    print(tab.pivot(index="condition", columns="method", values="detection_rate")
          .loc[ORDER, METHODS].round(2).to_string())
    print(f"\nPR anisotropy B = {pr_of(LAM_ANISO):.2f}; ring_s variance fraction in ring plane:",
          {s: round((1 + 2 * s * s) / (1 + D * s * s), 2) for s in SIGMAS},
          f"; K562 top-2 fraction {META['top2_var_frac']:.2f}")
    print(f"{time.time() - t:.0f}s")
