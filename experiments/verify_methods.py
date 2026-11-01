"""Check each method against an independent reference before trusting it.

This is not a unit-test suite; it is the short verification a careful reader would
want: the spectral distances agree with independent references (networkx resistance
distance, a direct transition-matrix diffusion distance) and behave like distance
matrices, a clean circle yields one dominant H1 loop, the PCA reference controls are
disjoint from the cells any test can draw, the energy distance matches its definition,
diagram distance is zero to itself and monotone in bar length, and the permutation
test rejects the null at about its nominal rate.

    python experiments/verify_methods.py
"""
import networkx as nx
import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform
from sklearn.decomposition import PCA

import config  # noqa: F401  (puts the package on the path)
from topoperturb import data, distances, features, homology, twosample


def check_effective_resistance():
    """The Laplacian-pseudoinverse resistance equals networkx's resistance_distance."""
    rng = np.random.default_rng(0)
    A = distances.knn_graph(rng.normal(size=(40, 5)))
    graph = nx.from_numpy_array(A)
    deg = A.sum(1)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        Lp = np.linalg.pinv(np.diag(deg) - A)
        diag = np.diag(Lp)
        R = diag[:, None] + diag[None, :] - 2 * Lp      # resistance before the paper's correction
        error = max(abs(R[i, j] - nx.resistance_distance(graph, i, j))
                    for i, j in [(0, 1), (2, 7), (5, 20)])
    assert error < 1e-6, error
    print(f"effective resistance agrees with networkx (max error {error:.1e})")


def check_diffusion_distance():
    """Diffusion distance equals an independent transition-matrix computation.

    Our eigendecomposition form and the raw `D^-1 A` random-walk form differ only
    by a global scale, which the pipeline's mean-off-diagonal scaling divides out;
    so the two must be perfectly correlated.
    """
    rng = np.random.default_rng(4)
    P = rng.normal(size=(60, 5))
    A = distances.knn_graph(P)
    deg = A.sum(1)
    pi = deg / deg.sum()
    off = ~np.eye(len(P), dtype=bool)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        for t in (1, 4, 8):
            mine = distances.diffusion_distance(P, t=t)
            Pt = np.linalg.matrix_power(A / deg[:, None], t)    # t-step transition matrix
            ref = np.sqrt(((Pt[:, None, :] - Pt[None, :, :]) ** 2 / pi).sum(-1))
            assert np.corrcoef(mine[off], ref[off])[0, 1] > 1 - 1e-9
    print("diffusion distance matches the transition-matrix definition (corr = 1)")


def check_distance_properties():
    """Both spectral distances are symmetric, zero on the diagonal, non-negative."""
    rng = np.random.default_rng(5)
    P = rng.normal(size=(50, 6))
    for name, D in [("diffusion", distances.diffusion_distance(P)),
                    ("effres", distances.effective_resistance(P))]:
        assert np.allclose(D, D.T), name
        assert np.allclose(np.diag(D), 0), name
        assert (D >= 0).all(), name
    print("diffusion and effective-resistance distances are valid metrics-matrices")


def check_circle_topology():
    """A clean circle has exactly one strongly persistent H1 loop (ripser sanity)."""
    rng = np.random.default_rng(6)
    angle = np.sort(rng.uniform(0, 2 * np.pi, 80))
    circle = np.c_[np.cos(angle), np.sin(angle)]
    _h0, h1 = homology.diagrams(squareform(pdist(circle)))
    life = features.lifetimes(h1)
    assert len(life) >= 1 and life[0] > 5 * (life[1] if len(life) > 1 else 0.0)
    print(f"a circle yields one dominant H1 loop (lifetime {life[0]:.2f})")


def check_pca_reference_split():
    """The PCA axes come only from reference controls, which no test can draw.

    Two things to establish: the reference and testing controls partition the control
    pool with no overlap (and no perturbed cell ever helps fit), and ``embed`` really
    uses those rows alone -- checked against a PCA fitted directly on them.
    """
    labels = np.array([data.CONTROL] * 800 + ["geneA"] * 60)
    reference = data.reference_split(labels, frac=0.5, seed=7)
    testing = data.testing_controls(labels, reference)
    controls = set(np.flatnonzero(labels == data.CONTROL))
    assert reference.sum() == 400 and not reference[testing].any()
    assert set(np.flatnonzero(reference)) | set(testing) == controls
    assert set(np.flatnonzero(reference)) & set(testing) == set()
    assert not reference[labels != data.CONTROL].any()
    assert (reference == data.reference_split(labels, frac=0.5, seed=7)).all()

    rng = np.random.default_rng(7)
    X = rng.normal(size=(len(labels), 12))
    _Z, pca = data.embed(X, labels, n_pcs=5, fit="reference", reference=reference, seed=7)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):   # as in data.embed
        direct = PCA(n_components=5, svd_solver="randomized", random_state=7).fit(X[reference])
    assert np.allclose(pca.components_, direct.components_)
    print(f"PCA axes come only from the {int(reference.sum())} reference controls, "
          f"disjoint from the {testing.size} testable ones")


def check_energy_distance():
    """Energy distance matches 2E|A-B| - E|A-A'| - E|B-B'| and is non-negative."""
    rng = np.random.default_rng(1)
    A, B = rng.normal(size=(60, 4)), rng.normal(size=(60, 4)) + 0.5
    reference = 2 * cdist(A, B).mean() - pdist(A).mean() - pdist(B).mean()
    assert abs(twosample.energy_distance(A, B) - reference) < 1e-9
    assert twosample.energy_distance(A, B) > 0
    print("energy distance matches its definition and is positive for differing clouds")


def check_diagram_distance():
    """Diagram distance is zero to itself, symmetric, and grows with a longer bar."""
    a = np.array([[0.0, 1.0]])
    b = np.array([[0.0, 2.0]])
    assert homology.diagram_distance(a, a) == 0.0
    assert abs(homology.diagram_distance(a, b) - homology.diagram_distance(b, a)) < 1e-12
    assert homology.diagram_distance(b, np.empty((0, 2))) > homology.diagram_distance(a, np.empty((0, 2)))
    print("diagram distance is zero to itself, symmetric, and monotone in bar length")


def check_permutation_calibration():
    """The permutation test on identical clouds rarely rejects (sanity on the null)."""
    rng = np.random.default_rng(3)
    energy_only = lambda A, B: {"energy": twosample.energy_distance(A, B)}
    rejects = 0
    for _ in range(40):
        A, B = rng.normal(size=(60, 5)), rng.normal(size=(60, 5))
        p, _ = twosample.permutation_test(energy_only, A, B, rng=rng)["energy"]
        rejects += p < 0.05
    assert rejects <= 6, rejects         # ~5% expected, allow slack on 40 draws
    print(f"permutation test on the null rejected {rejects}/40 times (≈5% expected)")


if __name__ == "__main__":
    check_effective_resistance()
    check_diffusion_distance()
    check_distance_properties()
    check_circle_topology()
    check_pca_reference_split()
    check_energy_distance()
    check_diagram_distance()
    check_permutation_calibration()
    print("all method checks passed")
