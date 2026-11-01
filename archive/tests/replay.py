"""Replay five archived comparisons, keeping the original 199/99 null draws.

Run from the repository root. Writes only results/local/replay.json; no downloads.
This is a spot check, not a new experiment or a full benchmark rerun.
"""
import importlib.metadata
import json
import os
import platform
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/topoperturb-matplotlib")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import synthetic as syn
import real
from methods import METHODS


def compare(actual, expected):
    fields = METHODS + [k for k in actual if k.startswith("obs_")]
    differences = {k: {"archived": float(expected[k]), "replayed": float(actual[k])}
                   for k in fields if not np.isclose(actual[k], expected[k], rtol=1e-5, atol=1e-7)}
    return {"matches_pvalues_and_statistics": not differences, "differences": differences,
            "replayed": {k: float(actual[k]) for k in fields}}


def replay_with_warnings(function, task):
    # Keep numerical warnings in the evidence without flooding the console.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = function(task)
    counts = Counter(f"{w.category.__name__}: {w.message}" for w in caught)
    return actual, dict(counts)


if __name__ == "__main__":
    start = time.monotonic()
    evidence = {"python": platform.python_version(), "platform": platform.platform(),
                "versions": {p: importlib.metadata.version(p) for p in
                             ["numpy", "scipy", "pandas", "scikit-learn", "ripser", "persim"]},
                "comparisons": []}
    out = ROOT / "results/local/replay.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    synthetic = pd.read_csv(syn.OUT / "synthetic_draws.csv", keep_default_na=False)
    real_draws = pd.read_csv(real.OUT / "real_draws.csv")
    for cond in ["null", "scaling", "ring_k562"]:
        actual, warning_counts = replay_with_warnings(syn.run, (list(syn.CONDITIONS).index(cond), cond, 0))
        expected = synthetic[(synthetic.condition == cond) & (synthetic.rep == 0)].iloc[0]
        result = {"case": f"synthetic/{cond}/rep0", "warning_counts": warning_counts, **compare(actual, expected)}
        evidence["comparisons"].append(result)
        print(result["case"], "MATCH" if result["matches_pvalues_and_statistics"] else result["differences"], flush=True)
        out.write_text(json.dumps(evidence, indent=2) + "\n")
    for cl, symbol in [("K562", "RPL21"), ("RPE1", "CPSF3")]:
        gi = next(i for i, (line, gene) in enumerate(real.GROUPS)
                  if line == cl and real.SYMS[cl].get(gene) == symbol)
        actual, warning_counts = replay_with_warnings(real.run, (gi, 0))
        expected = real_draws[(real_draws.cell_line == cl) & (real_draws.group == symbol) & (real_draws.draw == 0)].iloc[0]
        result = {"case": f"real/{cl}/{symbol}/draw0", "warning_counts": warning_counts, **compare(actual, expected)}
        evidence["comparisons"].append(result)
        print(result["case"], "MATCH" if result["matches_pvalues_and_statistics"] else result["differences"], flush=True)
        out.write_text(json.dumps(evidence, indent=2) + "\n")
    evidence["elapsed_seconds"] = time.monotonic() - start
    out.write_text(json.dumps(evidence, indent=2) + "\n")
    if any(not c["matches_pvalues_and_statistics"] for c in evidence["comparisons"]):
        raise SystemExit(1)
