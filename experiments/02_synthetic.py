"""Step 2a -- a synthetic benchmark with known truth.

We draw pairs of clouds where we control what differs: nothing (null), the mean
(translation), the scale, the anisotropy, the tails, the population mixture, or
the actual topology (a ring, a ring in noise, two clusters). The base Gaussian
uses the measured K562 control covariance, so sizes are realistic. For each
condition we report how often each method rejects at p < 0.05, with a Wilson
interval. A good topology method should fire on the ring/cluster rows and stay
near 5% on the geometry-only rows.

    python experiments/02_synthetic.py            # 40 replicates per condition
    python experiments/02_synthetic.py 10         # quicker

Note: the noisy ring and the mixture still have full-dimensional support, so a
pure distribution test is *allowed* to notice them; the truth labels describe the
latent shape, not a guarantee that only topology changed.
"""
import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

import config
from topoperturb import data, plotting, twosample

(_clouds, META) = data.load_clouds(config.cloud_path("k562"))
SPECTRUM = META["control_spectrum"]
LAM = np.array(SPECTRUM["eigenvalues"])            # K562 control variance per PC
SHIFT = SPECTRUM["median_centroid_shift"]          # a realistic translation size
D = len(LAM)
N, SEED = config.N_CELLS, config.SEED
SIGMAS = [0.1, 0.2, 0.3, 0.5]

# condition -> (family, what-changed); order fixes the figure rows and the seeds
CONDITIONS = {
    "null": ("geometry", "nothing"),
    "translation": ("geometry", "mean"),
    "scaling": ("geometry", "scale"),
    "anisotropy": ("geometry", "covariance shape"),
    "heavy_tails": ("geometry", "tails (non-Gaussian)"),
    "partial_kd": ("heterogeneity", "mixture"),
    "ring": ("topology", "loop (H1)"),
    "bimodal": ("topology", "two clusters (H0)"),
    **{f"ring_noise_{s}": ("topology", "loop (H1)") for s in SIGMAS},
}


def participation_ratio(values):
    """Effective number of dimensions carrying the variance."""
    return values.sum() ** 2 / (values ** 2).sum()


def sharpen(values, target_pr):
    """Reshape an eigenspectrum to a target participation ratio, keeping total variance."""
    power = brentq(lambda g: participation_ratio(values ** g) - target_pr, 1.0, 5.0)
    out = values ** power
    return out * values.sum() / out.sum()


LAM_ANISO = sharpen(LAM, 10.0)


def gaussian(rng, variances, n=N):
    """A centred Gaussian cloud with the given per-axis variances."""
    return rng.normal(size=(n, len(variances))) * np.sqrt(variances)


def make(condition, rng):
    """Return a (control, perturbed) pair for one benchmark condition."""
    if condition.startswith("ring_noise_"):
        sigma = float(condition.split("_")[-1])
        angle = rng.uniform(0, 2 * np.pi, N)
        perturbed = rng.normal(scale=sigma, size=(N, D))
        perturbed[:, 0] += np.cos(angle)
        perturbed[:, 1] += np.sin(angle)
        base = np.full(D, sigma ** 2)
        base[:2] += 0.5
        return gaussian(rng, base), perturbed
    control = gaussian(rng, LAM)
    if condition == "null":
        perturbed = gaussian(rng, LAM)
    elif condition == "translation":
        direction = rng.normal(size=D)
        perturbed = gaussian(rng, LAM) + SHIFT * direction / np.linalg.norm(direction)
    elif condition == "scaling":
        perturbed = 1.3 * gaussian(rng, LAM)
    elif condition == "anisotropy":
        perturbed = gaussian(rng, LAM_ANISO)
    elif condition == "heavy_tails":
        df = 5
        weight = rng.chisquare(df, N) / df
        perturbed = gaussian(rng, LAM) * np.sqrt((df - 2) / df) / np.sqrt(weight)[:, None]
    elif condition == "ring":
        angle = rng.uniform(0, 2 * np.pi, N)
        perturbed = gaussian(rng, LAM)
        perturbed[:, 0] = np.sqrt(2 * LAM[0]) * np.cos(angle)
        perturbed[:, 1] = np.sqrt(2 * LAM[1]) * np.sin(angle)
    elif condition == "bimodal":
        gap = 0.9 * np.sqrt(LAM[0])
        perturbed = gaussian(rng, LAM)
        perturbed[:, 0] = rng.choice([-gap, gap], N) + rng.normal(scale=np.sqrt(LAM[0] - gap ** 2), size=N)
    elif condition == "partial_kd":
        direction = rng.normal(size=D)
        perturbed = gaussian(rng, LAM)
        perturbed[: int(0.7 * N)] += SHIFT * direction / np.linalg.norm(direction)
    return control, perturbed


def run(task):
    """One replicate of one condition: generate the pair and run the battery."""
    index, condition, rep = task
    rng = np.random.default_rng([SEED, index, rep])
    control, perturbed = make(condition, rng)
    return {"condition": condition, "rep": rep, **twosample.compare(control, perturbed, rng)}


def summarise(draws):
    """Detection rate (p < 0.05) per condition and method, with Wilson intervals."""
    rows = []
    for condition, (family, changed) in CONDITIONS.items():
        subset = draws[draws.condition == condition]
        for method in twosample.METHODS:
            k = int((subset[f"p_{method}"] < 0.05).sum())
            lo, hi = twosample.wilson_interval(k, len(subset))
            rows.append({"condition": condition, "family": family, "changed": changed,
                         "method": method, "detection_rate": k / len(subset),
                         "ci_low": lo, "ci_high": hi, "n_reps": len(subset)})
    return pd.DataFrame(rows)


def main(repeats, out):
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(i, c, r) for i, c in enumerate(CONDITIONS) for r in range(repeats)]
    with Pool(8) as pool:
        draws = pd.DataFrame(pool.map(run, tasks, chunksize=1))
    draws.to_csv(out / "synthetic_draws.csv", index=False)
    table = summarise(draws)
    table.to_csv(out / "synthetic.csv", index=False)

    matrix = table.pivot(index="condition", columns="method", values="detection_rate")
    matrix = matrix.loc[list(CONDITIONS), twosample.METHODS]
    fig, ax = plotting.plt.subplots(figsize=(7.5, 6))
    plotting.heatmap(ax, matrix.values,
                     [f"{c}  ({CONDITIONS[c][1]})" for c in matrix.index],
                     twosample.METHODS,
                     "Synthetic benchmark: % of replicates with p < 0.05")
    plotting.save(fig, out / "synthetic.png")
    print(matrix.round(2).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repeats", type=int, nargs="?", default=40)
    parser.add_argument("--output", type=Path, default=config.RESULTS)
    args = parser.parse_args()
    t = time.time()
    main(args.repeats, args.output)
    print(f"{time.time() - t:.0f}s")
