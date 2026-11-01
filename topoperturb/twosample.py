"""Two-sample tests comparing a control cloud with a perturbed cloud.

Three questions, kept separate on purpose (see ``docs/METHODS.md``):

1. Do the two clouds differ at all?        -> permutation test of any statistic
2. Does the difference exceed a Gaussian?  -> covariance-matched surrogate test
3. What changed about the shape?           -> see ``features.py``

A statistic is any ``f(A, B) -> float`` that is larger when the clouds differ.
``permutation_test`` and ``surrogate_test`` work for a whole dict of statistics at
once, so every statistic is judged under the same random draws -- correct and
efficient. Persistent-homology statistics reuse one diagram summary per cloud.
"""
import numpy as np
from scipy.spatial.distance import pdist

from . import homology

P_PERM = 199      # label shuffles; p-value resolution is 1 / (P_PERM + 1)
M_SURR = 99       # Gaussian surrogate pairs

# Each persistent-homology feature is one (representation, dimension) pair.
PH_FEATURES = [("euclidean", "H0"), ("euclidean", "H1"),
               ("diffusion", "H1"), ("effres", "H1")]
PH_NAMES = [f"{rep}_{dim}" for rep, dim in PH_FEATURES]
METHODS = ["energy", "covariance"] + PH_NAMES


def energy_distance(A, B):
    """Energy distance between two clouds (Szekely & Rizzo); 0 iff identical."""
    cross = np.linalg.norm(A[:, None] - B[None], axis=2).mean()
    return 2 * cross - pdist(A).mean() - pdist(B).mean()


def cov_frobenius(A, B):
    """Frobenius norm of the covariance difference: a pure second-moment contrast."""
    return float(np.linalg.norm(np.cov(A.T) - np.cov(B.T)))


def ph_statistics(A, B):
    """Diagram distance per persistent-homology feature between two clouds."""
    reps = sorted({rep for rep, _ in PH_FEATURES})
    sA, sB = homology.summarize(A, reps), homology.summarize(B, reps)
    return {f"{rep}_{dim}": homology.diagram_distance(sA[rep][dim], sB[rep][dim])
            for rep, dim in PH_FEATURES}


def all_statistics(A, B):
    """The full method set: conventional contrasts plus the homology features."""
    return {"energy": energy_distance(A, B), "covariance": cov_frobenius(A, B),
            **ph_statistics(A, B)}


def wilson_interval(k, n, z=1.96):
    """95% Wilson confidence interval for a detection rate of ``k`` out of ``n``.

    Used to report how uncertain a detection or false-positive rate is.
    """
    if n == 0:
        return 0.0, 1.0
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def _pvalue(count, n):
    """Permutation p-value with the usual +1 so it is never exactly zero."""
    return (1 + count) / (1 + n)


def permutation_test(statistics, A, B, n=P_PERM, rng=None):
    """Shuffle the pooled cells ``n`` times; ask how often chance beats observed.

    ``statistics(X, Y)`` returns a ``{name: value}`` dict. Returns
    ``{name: (p_value, observed)}``. One shared set of shuffles tests every
    statistic, which also lets homology statistics reuse their diagrams.
    """
    rng = rng or np.random.default_rng()
    observed = statistics(A, B)
    pool, n_a = np.vstack([A, B]), len(A)
    exceed = dict.fromkeys(observed, 0)
    for _ in range(n):
        order = rng.permutation(len(pool))
        values = statistics(pool[order[:n_a]], pool[order[n_a:]])
        for name, value in values.items():
            exceed[name] += value >= observed[name]
    return {name: (_pvalue(exceed[name], n), observed[name]) for name in observed}


def surrogate_test(statistics, A, B, n=M_SURR, rng=None):
    """Compare the observed statistic with Gaussian clouds matched to each group.

    Each surrogate pair is drawn from Normals with the two groups' own means and
    covariances. Rejection means the difference is bigger than mean and covariance
    alone would produce -- it does *not* by itself prove a loop; any non-Gaussian
    structure can cause it. Returns ``{name: (p_value, observed, z_score)}``.
    """
    rng = rng or np.random.default_rng()
    observed = statistics(A, B)
    mean_a, cov_a = A.mean(0), np.cov(A.T)
    mean_b, cov_b = B.mean(0), np.cov(B.T)
    draws = {name: [] for name in observed}
    # A structured cloud (e.g. a ring) has a near-singular covariance; sampling
    # from it is valid but numerically noisy, so quiet the expected warnings.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for _ in range(n):
            A_ = rng.multivariate_normal(mean_a, cov_a, len(A))
            B_ = rng.multivariate_normal(mean_b, cov_b, len(B))
            for name, value in statistics(A_, B_).items():
                draws[name].append(value)
    results = {}
    for name, obs in observed.items():
        d = np.array(draws[name])
        p = (1 + (d >= obs).sum()) / (1 + n)
        spread = d.std(ddof=1)
        z = float((obs - d.mean()) / spread) if spread > 0 else 0.0  # all surrogates equal
        results[name] = (p, obs, z)
    return results


def compare(A, B, rng, n_perm=P_PERM, n_surr=M_SURR):
    """Run every test for one control/perturbed pair and flatten the results.

    Columns: ``p_<method>`` and ``obs_<method>`` for the permutation tests, and
    ``surr_<feature>`` / ``z_<feature>`` for the surrogate tests on each PH feature.
    """
    row = {}
    for name, (p, obs) in permutation_test(all_statistics, A, B, n_perm, rng).items():
        row[f"p_{name}"], row[f"obs_{name}"] = p, obs
    for name, (p, _obs, z) in surrogate_test(ph_statistics, A, B, n_surr, rng).items():
        row[f"surr_{name}"], row[f"z_{name}"] = p, z
    return row
