"""Recovered inputs must keep cell/gene order and be filtered consistently."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from prepare_data import load_cache


def test_cache_preserves_layout_and_filters_labels(tmp_path):
    path = tmp_path / "expression.npz"
    X = np.asfortranarray(np.arange(8, dtype=np.float32).reshape(4, 2))
    labels = np.array(["non-targeting", "gene-a", "gene-b", "excluded"])
    np.savez_compressed(path, X=X, lab=labels, genes=["gene-a", "gene-b"])
    filtered, lab, genes = load_cache(path, {"non-targeting", "gene-b"})
    np.testing.assert_array_equal(filtered, X[[0, 2]])
    assert list(lab) == ["non-targeting", "gene-b"]
    assert genes == ["gene-a", "gene-b"]

    np.savez_compressed(path, X=np.asfortranarray(X[:3]), lab=labels[:3], genes=genes)
    untouched, _, _ = load_cache(path)
    np.testing.assert_array_equal(untouched, X[:3])
    assert untouched.flags.f_contiguous


def test_cache_rejects_inconsistent_or_missing_inputs(tmp_path):
    path = tmp_path / "expression.npz"
    np.savez(path, X=np.ones((2, 2)), lab=["a"], genes=["a", "b"])
    with pytest.raises(ValueError, match="shapes"):
        load_cache(path)
    np.savez(path, X=np.ones((2, 2)), lab=["a", "b"], genes=["a", "b"])
    with pytest.raises(ValueError, match="missing requested"):
        load_cache(path, {"non-targeting"})
