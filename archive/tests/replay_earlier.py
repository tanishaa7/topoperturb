"""Replay one original comparison per recovered analysis, without full grids.

Uses original 2,000-permutation settings. Requires recovered expression inputs.
Writes a JSON report and exits unsuccessfully if a comparison differs.
"""
import argparse
from collections import Counter
import io
import json
import os
from pathlib import Path
import sys
import warnings
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/local/earlier_check")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "phase2").mkdir(exist_ok=True)
    os.environ["TOPOPERTURB_RESULTS"] = str(args.output.resolve())
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/topoperturb-matplotlib")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    sys.path.insert(0, str(ROOT / "src/earlier"))
    import numpy as np
    import pandas as pd

    results = []

    def table(member):
        with zipfile.ZipFile(ROOT / "docs/archive.zip") as z:
            return pd.read_csv(io.BytesIO(z.read(member)))

    def check(name, actual, expected, keys):
        differences = {k: {"replayed": float(actual[k]), "archived": float(expected[k])}
                       for k in keys if not np.isclose(actual[k], expected[k], rtol=1e-5, atol=1e-6)}
        results.append({"case": name, "matches": not differences, "differences": differences})
        print(name, "MATCH" if not differences else differences, flush=True)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        import permutations
        import normalized_ph
        import baselines
        import rpe1

        actual = permutations.run_gene(0)
        expected = table("topoperturb_phase1_outputs/perm2000_significance.csv")
        expected = expected[expected.gene_symbol == actual["gene_symbol"]].iloc[0]
        check("phase1/POLRMT", actual, expected,
              ["wass_H0", "wass_H1", "perm_H1_mean", "perm_H1_p95", "perm_p_H0", "perm_p_H1"])

        actual = normalized_ph.run(("primary", "POLRMT", 30, 150, 0, "euclid", 2000))
        expected = table("topoperturb_phase2_outputs/p2_ph_raw_results.csv")
        expected = expected[(expected.kind == "primary") & (expected.gene_symbol == "POLRMT")].iloc[0]
        check("phase2/normalized/POLRMT", actual, expected, ["obs", "perm_mean", "perm_p95", "p"])

        rows, surrogate = baselines.run_gene(0)
        expected = table("topoperturb_phase2_outputs/conventional_measures.csv")
        for actual in rows:
            row = expected[(expected.gene_symbol == "POLRMT") & (expected.measure == actual["measure"])].iloc[0]
            check("phase2/baseline/POLRMT/" + actual["measure"], actual, row, ["obs", "perm_mean", "p"])
        expected = table("topoperturb_phase2_outputs/gaussian_surrogates.csv")
        check("phase2/gaussian/POLRMT", surrogate, expected[expected.gene_symbol == "POLRMT"].iloc[0],
              ["surr_ratio", "surr_wins", "surr_obs_mean", "surr_null_mean"])

        actual = rpe1.run_gene(0)
        expected = table("topoperturb_phase2_outputs/rpe1_validation.csv")
        expected = expected[expected.label == "ZMAT2"].iloc[0]
        check("phase2/RPE1/ZMAT2", actual, expected,
              ["centroid_shift", "centroid_null_p95", "p_centroid", "raw_H1", "norm_H1", "energy",
               "p_raw_H1", "p_norm_H1", "p_energy"])

    report = {"comparisons": results, "warning_counts": dict(Counter(str(w.message) for w in caught))}
    (args.output / "replay.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if all(r["matches"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
