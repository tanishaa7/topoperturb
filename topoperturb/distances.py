"""Graph distances that expose a point cloud's shape better than raw Euclidean.

Persistent homology sees only the distance matrix it is given. On a kNN graph of
the cells, two intrinsic distances recover loops and clusters that Euclidean
distance blurs: the diffusion distance and the effective-resistance distance.
Both follow Damrich, Berens & Kobak, "Persistent homology for high-dimensional
data based on spectral methods" (NeurIPS 2024).
"""
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import kneighbors_graph

K_NN = 15        # neighbours in the graph
T_DIFF = 8       # diffusion time (how far random walks spread)


def knn_graph(P, k=K_NN):
    """Symmetric k-nearest-neighbour adjacency, connected into one component.

    If the graph is disconnected the spectral distances are undefined between
    components, so k is raised until a single component remains.
    """
    while True:
        A = kneighbors_graph(P, k, mode="connectivity", include_self=False).toarray()
        A = np.maximum(A, A.T)
        if connected_components(A, directed=False)[0] == 1:
            return A
        k += 5


def effective_resistance(P, k=K_NN):
    """Effective-resistance distance on the kNN graph (von Luxburg corrected).

    Resistance distance behaves like commute time between nodes; the correction
    removes the degree-driven term that otherwise dominates on large graphs.
    """
    A = knn_graph(P, k)
    deg = A.sum(1)
    # Near-degenerate clouds (e.g. Gaussian surrogates of a flat ring) make the
    # pseudoinverse numerically touchy; the clip below keeps the result valid.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        Lp = np.linalg.pinv(np.diag(deg) - A)        # pseudoinverse of the Laplacian
        diag = np.diag(Lp)
        R = diag[:, None] + diag[None, :] - 2 * Lp
        R = R - 1 / deg[:, None] - 1 / deg[None, :] + 2 * A / np.outer(deg, deg)
    R = np.clip(np.nan_to_num(R), 0, None)
    np.fill_diagonal(R, 0)
    return R


def diffusion_distance(P, k=K_NN, t=T_DIFF):
    """Diffusion distance on the kNN graph after ``t`` steps of a random walk.

    Cells are close when a random walk reaches similar places, which follows the
    manifold rather than straight-line gaps. Computed from the eigendecomposition
    of the symmetric normalised adjacency; the constant first eigenvector is dropped.
    """
    A = knn_graph(P, k)
    deg = A.sum(1)
    s = 1 / np.sqrt(deg)
    with np.errstate(over="ignore", invalid="ignore"):
        values, vectors = np.linalg.eigh(s[:, None] * A * s[None, :])
        order = np.argsort(values)[::-1]
        values, vectors = values[order], vectors[:, order]
        coords = (s[:, None] * vectors)[:, 1:] * (values[1:] ** t)
    return squareform(pdist(np.nan_to_num(coords)))
