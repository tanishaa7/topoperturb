"""Persistent homology of a point cloud, via ripser, and diagram comparison.

A persistence diagram records when topological features appear and disappear as a
ball of growing radius is drawn around the points: H0 points track connected
components (clusters), H1 points track loops. Each point is a (birth, death) pair;
its lifetime ``death - birth`` is how persistent -- and so how real -- the feature is.

The actual computation is delegated to established libraries: ``ripser`` builds the
Vietoris-Rips persistence, ``persim`` measures the distance between two diagrams.
"""
import numpy as np
from persim import wasserstein
from ripser import ripser
from scipy.spatial.distance import pdist, squareform

from . import distances

# Each representation gives its own diagram; a representation + "H1"/"H0" names a feature.
REPRESENTATIONS = ("euclidean", "diffusion", "effres")


def _finite(diagram):
    """Drop the one infinite H0 bar so diagrams only hold real (birth, death) pairs."""
    diagram = diagram[np.isfinite(diagram[:, 1])]
    return diagram if len(diagram) else np.empty((0, 2))


def diagrams(D):
    """Return the H0 and H1 persistence diagrams for a distance matrix ``D``."""
    both = ripser(D, maxdim=1, distance_matrix=True)["dgms"]
    return _finite(both[0]), _finite(both[1])


def _scaled(D):
    """Divide a distance matrix by its mean off-diagonal entry.

    This makes a cloud's diagram comparable across clouds of different overall
    scale, so a size difference alone does not look like a shape difference.
    """
    off = squareform(D, checks=False)      # the off-diagonal distances
    mean = off.mean() if off.size else 0.0
    return D / mean if mean > 1e-12 else D  # a near-degenerate cloud needs no rescaling


def representation(P, name):
    """Scaled distance matrix for one representation of the cloud ``P``."""
    if name == "euclidean":
        D = squareform(pdist(P))
    elif name == "diffusion":
        D = distances.diffusion_distance(P)
    elif name == "effres":
        D = distances.effective_resistance(P)
    else:
        raise ValueError(f"unknown representation {name!r}")
    return _scaled(D)


def diagram_distance(a, b):
    """Wasserstein distance between two persistence diagrams (0 if both empty)."""
    if len(a) == 0 and len(b) == 0:
        return 0.0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return float(wasserstein(a if len(a) else np.empty((0, 2)),
                                 b if len(b) else np.empty((0, 2))))


def summarize(P, reps=REPRESENTATIONS):
    """Diagrams for a cloud under each representation: ``{rep: {"H0", "H1"}}``.

    Computed once per cloud so several features can reuse the same diagrams.
    """
    out = {}
    for name in reps:
        h0, h1 = diagrams(representation(P, name))
        out[name] = {"H0": h0, "H1": h1}
    return out
