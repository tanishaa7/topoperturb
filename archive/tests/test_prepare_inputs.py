"""Scientific invariants for the new prepared-data route."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from prepare_inputs import audit, fit_pca, load_selected, select_targets
from fetch_inputs import fetch, SHA256


def test_loading_preserves_values_order_and_excludes_unsupported(tmp_path):
    # Cross a chunk boundary; row identity must remain global, not restart at 0.
    n = 10003
    labels = np.full(n, "excluded", dtype="U16")
    labels[[1, 10000]] = "non-targeting"
    labels[[3, 10002]] = "g1"
    table = pd.DataFrame({"g1": np.arange(n, dtype=np.float32) / 100,
                          "g2": np.ones(n), "intervention": labels})
    path = tmp_path / "input.tsv.gz"
    table.to_csv(path, sep="\t", index=False)
    genes, counts, report = audit(path)
    assert report["n_rows"] == n
    assert counts["excluded"] == n - 4
    X, lab, rows = load_selected(path, genes, ["g1"])
    np.testing.assert_array_equal(rows, [1, 3, 10000, 10002])
    np.testing.assert_allclose(X, table.loc[rows, genes].to_numpy(), rtol=1e-7)
    assert list(lab) == ["non-targeting", "g1", "non-targeting", "g1"]


def test_selection_uses_counts_and_never_excluded():
    from collections import Counter
    counts = Counter({"non-targeting": 400, "excluded": 900, "a": 150, "b": 149})
    assert select_targets(counts, None, True, 150, 300) == ["a"]
    for targets in (["excluded"], ["non-targeting"], ["missing"], ["b"]):
        with pytest.raises(ValueError):
            select_targets(counts, targets, False, 150, 300)


def test_balanced_pca_and_saved_transform_are_reproducible():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 6))
    X[:, -1] = 1
    labels = np.array(["non-targeting"] * 30 + ["a"] * 12 + ["b"] * 8)
    Z, pca, scale, fit, retained = fit_pca(X, labels, "balanced", 3, 10, "unit-variance", 42)
    assert dict(zip(*np.unique(labels[fit], return_counts=True))) == {"a": 8, "b": 8, "non-targeting": 8}
    assert not retained[-1]
    np.testing.assert_allclose(Z, ((X[:, retained] / scale) - pca.mean_) @ pca.components_.T)
    again = fit_pca(X, labels, "balanced", 3, 10, "unit-variance", 42)
    np.testing.assert_array_equal(Z, again[0])


def test_control_fit_does_not_use_perturbed_cells():
    X = np.random.default_rng(4).normal(size=(40, 5))
    labels = np.array(["non-targeting"] * 20 + ["a"] * 20)
    first = fit_pca(X, labels, "controls", 3, 10, "none", 42)
    X[20:] += 500
    second = fit_pca(X, labels, "controls", 3, 10, "none", 42)
    np.testing.assert_array_equal(first[1].components_, second[1].components_)
    np.testing.assert_allclose(first[0][:20], second[0][:20])
    with pytest.raises(ValueError, match="variable"):
        fit_pca(np.zeros((40, 5)), labels, "pooled", 3, 10, "none", 42)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1])
def test_audit_rejects_invalid_expression(tmp_path, bad):
    path = tmp_path / "bad.tsv.gz"
    pd.DataFrame({"g1": [bad], "intervention": ["non-targeting"]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="finite"):
        audit(path)


def test_fetch_refuses_unverified_existing_file(tmp_path):
    path = tmp_path / "causalbench_rpe1.tsv.gz"
    path.write_bytes(b"not the dataset")
    with pytest.raises(ValueError, match="pinned"):
        fetch("rpe1", tmp_path)
    assert path.read_bytes() == b"not the dataset"
