"""
Supplementary figure: distribution of eigenvalues for forest (real data) and
SENM (neutral model), at two scales, with the MP cutoff marked as a vertical
dotted line.

Run: python supp/distributions.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from src.compute import compute_forest_spectrum, compute_senm_spectrum, get_dimensions


OUT_DIR = "figures_supp"

LABELS = {'forest': 'BCI', 'senm': 'SENM'}
COLORS = {'forest': '#2ca02c', 'senm': '#1f77b4'}
SENM_XMAX = {7: 8, 45: 22}   # scale -> right limit of the SENM panels


def mp_cutoff(Q: float, sigma2: float = 1.0) -> float | None:
    """
    Upper MP edge lambda_+, in terms of Q = T/N (validity Q >= 1).
    Returns None outside the valid range.
    """
    if Q < 1:
        return None
    return sigma2 * (1 + np.sqrt(1 / Q)) ** 2


def compute_spectra_by_scale(
    scales: list[int],
    n_species: int,
    censuses: list[int],
    reference_census: int,
    n_censuses_senm: int,
) -> dict:
    """Compute forest and senm spectra once per scale, pooled across units."""
    spectra_by_scale = {}
    for scale in scales:
        print(f"Scale {scale}: computing forest and senm spectra...")
        _, _, forest_spectra = compute_forest_spectrum(
            n_species, scale, censuses, reference_census
        )
        _, _, senm_spectra, _ = compute_senm_spectrum(
            n_species, scale, n_censuses_senm
        )
        spectra_by_scale[scale] = {
            'forest': forest_spectra.ravel(),
            'senm': senm_spectra.ravel(),
        }
    return spectra_by_scale


def plot_eig_distributions(
    scales: list[int],
    n_species: int,
    censuses: list[int],
    reference_census: int,
    n_censuses_senm: int,
):
    os.makedirs(OUT_DIR, exist_ok=True)
    spectra_by_scale = compute_spectra_by_scale(
        scales, n_species, censuses, reference_census, n_censuses_senm
    )
    

    for data_type in ('forest', 'senm'):
        fig, axes = plt.subplots(1, len(scales), figsize=(6 * len(scales), 5))
        if len(scales) == 1:
            axes = [axes]

        grid_width, grid_height = get_dimensions(data_type)

        for ax, scale in zip(axes, scales):
            eigs = spectra_by_scale[scale][data_type]

            T = int(grid_width / scale) * int(grid_height / scale)
            Q = T / n_species
            lambda_max = mp_cutoff(Q)

            ax.hist(eigs, bins=40, density=True, alpha=0.6,
                    color=COLORS[data_type], label=LABELS[data_type])

            if lambda_max is not None:
                ax.axvline(lambda_max, color='crimson', ls=':', lw=2, label='MP cutoff')
            else:
                ax.text(0.5, 0.5, 'Q < 1: outside MP range', transform=ax.transAxes,
                        ha='center', va='center', color='crimson')

            ax.axvline(eigs.max(), color='black', ls='--', lw=1.2,
                       label=rf'$\lambda_{{\max}}$ = {eigs.max():.2f}')

            ax.set_title(f'Scale = {scale} m  (Q = {Q:.2f})')
            ax.set_xlabel('Eigenvalue')
            ax.set_ylabel('Density')
            ax.spines[['top', 'right']].set_visible(False)
            ax.legend(frameon=False, fontsize=9)
            if data_type == 'senm':
                ax.set_xlim(right=SENM_XMAX[scale])

        fig.suptitle(f'{LABELS[data_type]} eigenvalue distribution')
        fig.tight_layout()

        out_path = f"{OUT_DIR}/fig_eig_distributions_{data_type}.svg"
        fig.savefig(out_path)
        print(f"Saved {out_path}")


def main():
    n_species = 100
    reference_census = 4
    censuses = list(range(1, 9))
    n_censuses_senm = 100
    scales = [7, 45]

    plot_eig_distributions(scales, n_species, censuses, reference_census, n_censuses_senm)


if __name__ == "__main__":
    main()
