"""Small offline checks of mathematics and archived benchmark bookkeeping.

These do not establish biological validity or rerun the full benchmark.
"""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist, squareform
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "src"
sys.path.insert(0, str(BENCH))
import methods as tb

import synthetic as syn
import real

def test_rips_square_and_empty_diagrams():
    points = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    h0, h1 = tb._dgms(squareform(pdist(points)))
    np.testing.assert_allclose(h0, [[0, 1]] * 3)
    np.testing.assert_allclose(h1, [[1, np.sqrt(2)]], rtol=1e-6)
    empty = np.empty((0, 2))
    assert tb.wass(empty, empty) == 0
    assert tb.wass(empty, h1) == pytest.approx((np.sqrt(2) - 1) / np.sqrt(2), abs=1e-6)


def test_scale_and_translation_invariance_of_statistics():
    p = np.random.default_rng(12).normal(size=(45, 6))
    a, b = tb.summary(p), tb.summary(1.3 * p + 5)
    np.testing.assert_allclose(b["raw_H1"], 1.3 * a["raw_H1"], rtol=2e-6)
    for k in tb.PH_KEYS[1:]:
        assert tb.wass(a[k], b[k]) < 1e-5


def test_diffusion_against_direct_random_walk():
    p = np.random.default_rng(2).normal(size=(45, 5))
    adj = tb.knn_adj(p)
    degree = adj.sum(1)
    walk = np.linalg.matrix_power(adj / degree[:, None], tb.T_DIFF)
    direct = squareform(pdist(walk / np.sqrt(degree)[None, :]))
    # The conventional sqrt(graph volume) factor cancels on mean normalization.
    np.testing.assert_allclose(tb.dist_diffusion(p), direct, rtol=1e-9, atol=1e-12)


def test_corrected_resistance_against_spectral_formula():
    p = np.random.default_rng(3).normal(size=(45, 5))
    adj = tb.knn_adj(p)
    degree = adj.sum(1)
    s = 1 / np.sqrt(degree)
    mu, u = np.linalg.eigh(np.eye(len(p)) - s[:, None] * adj * s[None, :])
    coords = s[:, None] * u[:, 1:] * ((1 - mu[1:]) / np.sqrt(mu[1:]))
    direct = squareform(pdist(coords, metric="sqeuclidean"))
    np.testing.assert_allclose(tb.dist_effres(p), direct, rtol=1e-8, atol=1e-12)


def test_baselines_against_explicit_u_statistic():
    rng = np.random.default_rng(9)
    a, b = rng.normal(size=(20, 4)), rng.normal(size=(25, 4))
    da, db = squareform(pdist(a)), squareform(pdist(b))
    cross = np.linalg.norm(a[:, None] - b[None, :], axis=2)
    expected = 2 * cross.mean() - da.sum() / (20 * 19) - db.sum() / (25 * 24)
    assert tb.energy(a, b) == pytest.approx(expected)
    ca = (a - a.mean(0)).T @ (a - a.mean(0)) / 19
    cb = (b - b.mean(0)).T @ (b - b.mean(0)) / 24
    assert tb.cov_frob(a, b) == pytest.approx(np.sqrt(((ca - cb) ** 2).sum()))


def test_permutation_plus_one_and_ties():
    # All statistics tie. With >= and a plus-one correction every p must be one.
    p = np.random.default_rng(3).normal(size=(20, 3))
    empty = {k: np.empty((0, 2)) for k in tb.PH_KEYS}
    with patch.object(tb, "summary", return_value=empty), \
         patch.object(tb, "energy", return_value=0.0), \
         patch.object(tb, "cov_frob", return_value=0.0):
        out = tb.battery(p, p, np.random.default_rng(4), n_perm=9, n_surr=3)
    assert all(out[k] == 1 for k in tb.METHODS)


def test_bh_against_scipy():
    p = np.array([0.5, 0.005, 0.005, 0.04, 1.0, 0.03])
    np.testing.assert_allclose(real.bh(p), false_discovery_control(p))


def test_synthetic_construction():
    np.testing.assert_allclose(syn.LAM_ANISO.sum(), syn.LAM.sum())
    assert syn.pr_of(syn.LAM_ANISO) == pytest.approx(10)
    for cond in syn.CONDITIONS:
        a, b = syn.make(cond, np.random.default_rng(18))
        assert a.shape == b.shape == (150, 30)
        assert np.isfinite(a).all() and np.isfinite(b).all()
    _, ring = syn.make("ring_k562", np.random.default_rng(18))
    np.testing.assert_allclose((ring[:, :2] ** 2 / (2 * syn.LAM[:2])).sum(1), 1)


def test_all_synthetic_aggregates():
    draws = pd.read_csv(syn.OUT / "synthetic_draws.csv", keep_default_na=False)
    assert len(draws) == 480
    assert not draws.duplicated(["condition", "rep"]).any()
    assert draws.groupby("condition").rep.apply(lambda x: set(x) == set(range(40))).all()
    expected = pd.read_csv(syn.OUT / "synthetic.csv", keep_default_na=False)
    actual = syn.summarise(draws)
    pd.testing.assert_frame_equal(actual, expected, check_exact=False, rtol=1e-12, atol=1e-12)


def test_all_real_aggregates():
    draws = pd.read_csv(real.OUT / "real_draws.csv")
    table = pd.read_csv(real.OUT / "real.csv")
    assert len(draws) == 225 and len(table) == 150
    assert not draws.duplicated(["cell_line", "group", "draw"]).any()
    assert draws.groupby(["cell_line", "group"]).draw.apply(lambda x: set(x) == set(range(15))).all()
    for row in table.itertuples():
        d = draws[(draws.cell_line == row.cell_line) & (draws.group == row.group)]
        assert row.detect_frac == pytest.approx((d[row.method] < 0.05).mean())
        assert row.median_p == pytest.approx(d[row.method].median())
        assert row.p_primary == pytest.approx(d.loc[d.draw == 0, row.method].iloc[0])
    genes = ~table.group.str.startswith(("control", "pos."))
    for _, t in table[genes].groupby(["cell_line", "method"]):
        np.testing.assert_allclose(t.q_primary, false_discovery_control(t.p_primary))
    assert table.loc[~genes, "q_primary"].isna().all()


def test_saved_pvalue_resolution():
    for name in ["synthetic_draws.csv", "real_draws.csv"]:
        d = pd.read_csv(syn.OUT / name, keep_default_na=False)
        for m in tb.METHODS:
            denom = 100 if m.startswith("surr_") else 200
            assert d[m].between(1 / denom, 1).all()
            np.testing.assert_allclose(d[m] * denom, np.round(d[m] * denom), atol=1e-10)


def test_cloud_counts_and_metadata():
    for cl, (z, lab) in real.CLOUDS.items():
        assert z.shape[1] == 30 and np.isfinite(z).all()
        inventory = pd.read_csv(ROOT / "data/cell_counts.csv")
        counts = inventory[inventory.cell_line == cl].set_index("intervention").n_cells
        assert set(lab) == set(real.SYMS[cl]) | {"non-targeting"}
        for label, n in pd.Series(lab).value_counts().items():
            assert n == counts[label] and n >= 150
    z, lab = real.CLOUDS["K562"]
    lam = np.linalg.eigvalsh(np.cov(z[lab == "non-targeting"].T))[::-1]
    np.testing.assert_allclose(lam, syn.LAM, rtol=1e-10)
    assert syn.pr_of(lam) == pytest.approx(syn.META["PR_ctrl"])


@pytest.mark.xfail(strict=True, reason="Existing kNN implementation does not cap k at n-1")
def test_small_cloud_graph_boundary():
    assert tb.knn_adj(np.random.default_rng(1).normal(size=(10, 3))).shape == (10, 10)


@pytest.mark.xfail(strict=True, reason="Existing scale normalization divides by zero for collapsed clouds")
def test_collapsed_cloud_rejected_clearly():
    with pytest.raises(ValueError, match="scale|distinct|degenerate"):
        tb.summary(np.zeros((20, 3)))
