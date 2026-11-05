"""Eigenvalue spectra (BCI and SENM) and the BCI correlation matrix at every scale.

Writes data/processed/scale_<l>/:
    forest_spectra.npy                                  (n_censuses, n_species), descending
    senm_spectra.npy                                    (n_realizations, n_species), descending
    forest_correlation_mean.npy, forest_correlation_std.npy
"""
import os
import numpy as np

from src.config import load_config
from src.compute import compute_forest_spectrum, compute_forest_C, compute_senm_spectrum

config     = load_config()
n_species  = config['n_species']
censuses   = config['censuses']
ref_census = config['reference_census']
n_real     = config['senm']['num_realizations']

for scale in range(5, 50):
    out = f"data/processed/scale_{scale}"
    os.makedirs(out, exist_ok=True)
    print(f"scale {scale} m", flush=True)

    *_, forest_spectra = compute_forest_spectrum(n_species, scale, censuses, ref_census)
    C_mean, C_std, *_  = compute_forest_C(n_species, scale, censuses, ref_census)
    _, _, senm_spectra, _ = compute_senm_spectrum(n_species, scale, n_real)

    np.save(f"{out}/forest_spectra.npy", forest_spectra)
    np.save(f"{out}/senm_spectra.npy", senm_spectra)
    np.save(f"{out}/forest_correlation_mean.npy", C_mean)
    np.save(f"{out}/forest_correlation_std.npy", C_std)
