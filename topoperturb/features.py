"""Interpretable topology features: not just *that* clouds differ, but *how*.

A permutation p-value says two clouds differ; a surrogate p-value says the
difference beats a Gaussian. Neither says what changed or in which direction.
These helpers read that off the persistence diagrams directly: how long-lived the
main loop is, whether it grew or shrank from control to perturbed, and -- by
bootstrapping -- how stable that reading is.
"""
import numpy as np

from . import homology


def lifetimes(diagram):
    """Feature lifetimes ``death - birth``, longest first (empty -> empty)."""
    if len(diagram) == 0:
        return np.empty(0)
    return np.sort(diagram[:, 1] - diagram[:, 0])[::-1]


def h1_summary(P, representation="diffusion"):
    """Scalar loop summary for a cloud: leading/second lifetime, count, total.

    ``leading_H1`` is the most persistent loop -- the headline "is there a ring".
    """
    h1 = homology.summarize(P, (representation,))[representation]["H1"]
    life = lifetimes(h1)
    return {"leading_H1": float(life[0]) if len(life) else 0.0,
            "second_H1": float(life[1]) if len(life) > 1 else 0.0,
            "n_loops": int(len(life)),
            "total_H1": float(life.sum())}


def signed_change(control, perturbed, representation="diffusion"):
    """Change in leading H1 lifetime from control to perturbed, with its sign.

    Positive means a loop became more persistent under perturbation (appeared or
    strengthened); negative means it weakened or vanished. A distance alone cannot
    tell these apart, so this is what makes a "loop appeared/disappeared" claim.
    """
    before = h1_summary(control, representation)["leading_H1"]
    after = h1_summary(perturbed, representation)["leading_H1"]
    return {"leading_H1_control": before, "leading_H1_perturbed": after,
            "leading_H1_change": after - before}


def bootstrap_leading_h1(P, representation="diffusion", n=200, rng=None):
    """Resample cells to put an uncertainty band on the leading H1 lifetime.

    Returns the mean and a 95% percentile interval, so a reported loop never comes
    without a sense of how reproducible it is across resamples of the same cells.
    """
    rng = rng or np.random.default_rng()
    values = np.array([
        h1_summary(P[rng.integers(0, len(P), len(P))], representation)["leading_H1"]
        for _ in range(n)])
    lo, hi = np.percentile(values, [2.5, 97.5])
    return {"leading_H1_mean": float(values.mean()),
            "leading_H1_lo": float(lo), "leading_H1_hi": float(hi)}
