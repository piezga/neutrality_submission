"""Matrix data shared by Fig. 1c and the filtered-matrix SI figure, read from data/processed."""
import numpy as np
from src.config import load_config
from src.compute import marchenko_pastur_bounds, marchenko_pastur_filter
from src.communities import L_detect_communities
from src.utils import load_C, top_species, get_dimensions


def forest_matrix(scale=20):
    cfg = load_config(); N = cfg['n_species']
    C = load_C('forest', scale)[0]
    names = sorted(top_species('forest', N, cfg['reference_census']))           # rows of C are alphabetical
    lam_plus = marchenko_pastur_bounds('forest', N, scale)[0]
    CMP = marchenko_pastur_filter(C, scale, 'forest')                            # keeps every mode above lambda_+
    w = np.linalg.eigvalsh(C); keep = w > lam_plus
    _, labels, order = L_detect_communities(CMP, cfg['communities']['n_communities'])
    W, H = get_dimensions('forest')
    return dict(C=C, CMP=CMP, names=names, lam_plus=lam_plus, w=w, keep=keep, K=int(keep.sum()),
                eig=np.sort(w)[::-1], order=order, labels=labels,
                bounds=np.where(np.diff(labels[order]))[0] + 0.5, T=(W // scale) * (H // scale))
