"""Test battery shared by the synthetic benchmark and the real-data check.

For two point clouds A and B (cells x PCs) it returns one p-value per method:

  permutation tests (label shuffling of the pooled cells, P_PERM shuffles)
    energy      energy distance between the two clouds
    cov_frob    Frobenius norm of the covariance difference
    raw_H1      H1 Wasserstein distance, Euclidean Rips, no rescaling
    norm_H1     same, each cloud divided by its own mean pairwise distance
    norm_H0     H0 Wasserstein distance (finite bars), rescaled as norm_H1
    effres_H1   H1 on corrected effective resistance (kNN graph), rescaled
    diff_H1     H1 on diffusion distance (kNN graph, t = 8), rescaled

  covariance-matched Gaussian surrogate tests (M_SURR surrogate pairs)
    surr_norm_H1, surr_norm_H0, surr_diff_H1
    The observed PH distance is compared with the PH distance between
    Gaussian clouds that have each sample's own mean and covariance.
    Rejection means: the PH difference is larger than mean + covariance
    alone would produce. It does NOT by itself mean a loop exists; any
    non-Gaussian structure (tails, clusters, loops) can cause it.

Spectral distances follow Damrich, Berens & Kobak (NeurIPS 2024).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import kneighbors_graph
from ripser import ripser
from persim import wasserstein

K_NN, T_DIFF = 15, 8
P_PERM, M_SURR = 199, 99
PH_KEYS = ["raw_H1", "norm_H1", "norm_H0", "effres_H1", "diff_H1"]
SURR_KEYS = ["norm_H1", "norm_H0", "diff_H1"]
METHODS = (["energy", "cov_frob"] + PH_KEYS + ["surr_" + k for k in SURR_KEYS])
LABELS = {"energy": "Energy distance", "cov_frob": "Covariance diff.",
          "raw_H1": "Raw H1", "norm_H1": "Scale-norm. H1", "norm_H0": "Scale-norm. H0",
          "effres_H1": "Eff.-resistance H1", "diff_H1": "Diffusion H1",
          "surr_norm_H1": "Surrogate-adj. H1", "surr_norm_H0": "Surrogate-adj. H0",
          "surr_diff_H1": "Surrogate-adj. diffusion H1"}


def knn_adj(P):
    k = K_NN
    while True:
        A = kneighbors_graph(P, k, mode="connectivity", include_self=False).toarray()
        A = np.maximum(A, A.T)
        if connected_components(A, directed=False)[0] == 1:
            return A
        k += 5


def dist_effres(P):
    A = knn_adj(P)
    deg = A.sum(1)
    Lp = np.linalg.pinv(np.diag(deg) - A)
    dg = np.diag(Lp)
    R = dg[:, None] + dg[None, :] - 2 * Lp
    R = R - 1 / deg[:, None] - 1 / deg[None, :] + 2 * A / np.outer(deg, deg)   # von Luxburg correction
    R = np.clip(R, 0, None)
    np.fill_diagonal(R, 0)
    return R


def dist_diffusion(P):
    A = knn_adj(P)
    deg = A.sum(1)
    s = 1 / np.sqrt(deg)
    lam, V = np.linalg.eigh(s[:, None] * A * s[None, :])
    order = np.argsort(lam)[::-1]
    lam, V = lam[order], V[:, order]
    coords = (s[:, None] * V)[:, 1:] * (lam[1:] ** T_DIFF)
    return squareform(pdist(coords))


def _dgms(D):
    d = ripser(D, maxdim=1, distance_matrix=True)["dgms"]
    h0 = d[0][np.isfinite(d[0][:, 1])]
    h1 = d[1] if len(d[1]) else np.empty((0, 2))
    return h0, h1


def _mean_off(D):
    return D[np.triu_indices(len(D), 1)].mean()


def summary(P, keys=PH_KEYS):
    out = {}
    De = squareform(pdist(P))
    s = _mean_off(De)
    h0, h1 = _dgms(De)
    out["raw_H1"], out["norm_H1"], out["norm_H0"] = h1, h1 / s, h0 / s
    if "diff_H1" in keys:
        Dd = dist_diffusion(P)
        out["diff_H1"] = _dgms(Dd / _mean_off(Dd))[1]
    if "effres_H1" in keys:
        Dr = dist_effres(P)
        out["effres_H1"] = _dgms(Dr / _mean_off(Dr))[1]
    return out


def wass(a, b):
    if len(a) == 0 and len(b) == 0:
        return 0.0
    return float(wasserstein(a if len(a) else np.empty((0, 2)), b if len(b) else np.empty((0, 2))))


def energy(A, B):
    cross = np.linalg.norm(A[:, None] - B[None], axis=2).mean()
    return 2 * cross - pdist(A).mean() - pdist(B).mean()


def cov_frob(A, B):
    return float(np.linalg.norm(np.cov(A.T) - np.cov(B.T)))

def _stats(A, B, sa, sb, keys):
    out = {k: wass(sa[k], sb[k]) for k in keys}
    return out


def battery(A, B, rng, n_perm=P_PERM, n_surr=M_SURR):
    sa, sb = summary(A), summary(B)
    obs = _stats(A, B, sa, sb, PH_KEYS)
    obs["energy"], obs["cov_frob"] = energy(A, B), cov_frob(A, B)

    pool, na = np.vstack([A, B]), len(A)
    ge = dict.fromkeys(obs, 0)
    for _ in range(n_perm):
        q = rng.permutation(len(pool))
        X1, X2 = pool[q[:na]], pool[q[na:]]
        st = _stats(X1, X2, summary(X1), summary(X2), PH_KEYS)
        st["energy"], st["cov_frob"] = energy(X1, X2), cov_frob(X1, X2)
        for k, v in st.items():
            ge[k] += v >= obs[k]
    res = {k: (1 + ge[k]) / (1 + n_perm) for k in obs}

    ma, Ca, mb, Cb = A.mean(0), np.cov(A.T), B.mean(0), np.cov(B.T)
    sur = {k: [] for k in SURR_KEYS}
    for _ in range(n_surr):
        A_ = rng.multivariate_normal(ma, Ca, len(A))
        B_ = rng.multivariate_normal(mb, Cb, len(B))
        st = _stats(A_, B_, summary(A_, SURR_KEYS), summary(B_, SURR_KEYS), SURR_KEYS)
        for k in SURR_KEYS:
            sur[k].append(st[k])
    for k in SURR_KEYS:
        s = np.array(sur[k])
        res["surr_" + k] = (1 + (s >= obs[k]).sum()) / (1 + n_surr)
        res["z_surr_" + k] = (obs[k] - s.mean()) / s.std(ddof=1)
    for k in PH_KEYS + ["energy", "cov_frob"]:
        res["obs_" + k] = obs[k]
    return res
