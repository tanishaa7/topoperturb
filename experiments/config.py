"""Shared settings for the experiments: paths, the cells we study, and names.

The perturbations are the eight K562 pilot knockdowns, and for RPE1 the four
knockdowns shared with K562 plus one strong positive control (the largest
centroid shift among eligible genes). Genes are Ensembl IDs; symbols are for
display only.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))                 # import the topoperturb package when run as a script

DATA = ROOT / "data"
RESULTS = ROOT / "results"
N_CELLS, N_PCS, SEED = 150, 30, 42

TARGETS = {
    "k562": ["ENSG00000099821", "ENSG00000168827", "ENSG00000146007", "ENSG00000119203",
             "ENSG00000122026", "ENSG00000047315", "ENSG00000091140", "ENSG00000114503"],
    "rpe1": ["ENSG00000146007", "ENSG00000119203", "ENSG00000099821", "ENSG00000168827",
             "ENSG00000174243"],
}
SYMBOLS = {
    "ENSG00000099821": "POLRMT", "ENSG00000168827": "GFM1", "ENSG00000146007": "ZMAT2",
    "ENSG00000119203": "CPSF3", "ENSG00000122026": "RPL21", "ENSG00000047315": "POLR2B",
    "ENSG00000091140": "DLD", "ENSG00000114503": "NCBP2", "ENSG00000174243": "pos. control",
}


def symbol(gene):
    """Readable label for a gene ID (falls back to the ID itself)."""
    return SYMBOLS.get(gene, gene)


def data_source(cache):
    """Describe where the expression came from, for the saved provenance record."""
    return f"local cache {cache}" if cache else "Hugging Face (pinned revision)"


def cloud_path(cell_line, n_pcs=None):
    """Where ``01_prepare_data`` writes the PCA clouds for a cell line.

    The primary ``N_PCS`` representation keeps the plain name; other PC counts are
    tagged, so a representation-sensitivity run never overwrites the primary clouds.
    """
    name = "clouds.npz" if n_pcs in (None, N_PCS) else f"clouds_{n_pcs}pc.npz"
    return DATA / cell_line / name
