"""TopoPerturb: persistent homology of single-cell CRISPR perturbations.

A small, readable pipeline that asks whether knocking down a gene changes the
*shape* of the cell population -- a loop or cluster structure -- beyond what a
shift in mean or covariance explains. The parts:

    data        load CausalBench from Hugging Face and build PCA clouds
    distances   diffusion and effective-resistance distances on a kNN graph
    homology    persistence diagrams (ripser) and diagram distances (persim)
    twosample   energy / permutation / Gaussian-surrogate tests
    features    interpretable loop summaries with direction and uncertainty
    plotting    figures

The experiments in ``experiments/`` compose these into the synthetic benchmark,
the real knockdown analysis and the negative-control calibration.
"""
from . import data, distances, features, homology, plotting, twosample

__all__ = ["data", "distances", "features", "homology", "plotting", "twosample"]
