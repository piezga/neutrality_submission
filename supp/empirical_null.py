"""
Compare the empirical "spin" null model against the theoretical
Marchenko-Pastur distribution, at a couple of spatial scales.

For each scale:
  1. Build the real BCI spatial_series (n_species x n_bins) at that scale.
  2. Repeat 100 times: circularly shift each row by an independent random
     amount ("spin"), compute the correlation matrix, get its eigenvalues.
  3. Pool all eigenvalues across the 100 spins into one empirical distribution.
  4. Overlay the theoretical MP density (same N, T -> same Q) on top.

Run: python supp/empirical_null.py
"""

import numpy as np
import matplotlib.pyplot as plt

from src.utils import raw2df, top_species
from src.compute import compute_spatial_series, get_dimensions


# --------------------------------------------------------------------------
# Spin null model
# --------------------------------------------------------------------------

def spin_rows(spatial_series: np.ndarray, rng: np.random.Generator | None = None) -> np.ndarray:
    """
    Circularly shift each row of spatial_series by an independent random
    amount ("bike lock" spin), preserving each row's internal (spatial)
    structure while destroying cross-row correlations.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_species, n_bins = spatial_series.shape
    shifts = rng.integers(0, n_bins, size=n_species)

    idx = (np.arange(n_bins) - shifts[:, None]) % n_bins
    return np.take_along_axis(spatial_series, idx, axis=1)


def spin_columns(spatial_series: np.ndarray, rng: np.random.Generator | None = None) -> np.ndarray:
    """
    Circularly shift each column of spatial_series by an independent random
    amount along the species axis ("vertical" spin), destroying each
    species' own internal (spatial) structure while permuting which species
    identity sits at which rank, for every bin independently.

    Parameters
    ----------
    spatial_series : np.ndarray
        Shape (n_species, n_bins).
    rng : np.random.Generator, optional

    Returns
    -------
    np.ndarray, shape (n_species, n_bins)
    """
    if rng is None:
        rng = np.random.default_rng()

    n_species, n_bins = spatial_series.shape
    shifts = rng.integers(0, n_species, size=n_bins)

    idx = (np.arange(n_species)[:, None] - shifts[None, :]) % n_species
    return np.take_along_axis(spatial_series, idx, axis=0)

def spin_null_eigenvalues(
    spatial_series: np.ndarray,
    n_spins: int = 100,
    rng: np.random.Generator | None = None,
    mode: str = 'rows',
) -> np.ndarray:
    """
    Pool eigenvalues from n_spins independent spins of spatial_series.

    mode : str
        'rows'    -- row-spin only: preserves each species' internal
                     spatial (autocorrelation) structure, destroys
                     cross-species correlation.
        'columns' -- column-spin only: destroys each species' internal
                     spatial structure, isolating what happens to
                     cross-species correlation alone.
        'both'    -- row-spin AND column-spin combined: destroys both,
                     as a maximally-scrambled baseline.

    Returns
    -------
    np.ndarray, shape (n_spins * n_species,)
    """
    if mode not in ('rows', 'columns', 'both'):
        raise ValueError(f"mode must be 'rows', 'columns', or 'both', got {mode!r}")

    if rng is None:
        rng = np.random.default_rng()

    n_species, n_bins = spatial_series.shape
    all_eigs = np.empty((n_spins, n_species))

    for i in range(n_spins):
        spun = spatial_series
        if mode in ('rows', 'both'):
            spun = spin_rows(spun, rng=rng)
        if mode in ('columns', 'both'):
            spun = spin_columns(spun, rng=rng)
        C = np.corrcoef(spun)
        all_eigs[i] = np.linalg.eigvalsh(C)

    return all_eigs.ravel()
# --------------------------------------------------------------------------
# Theoretical Marchenko-Pastur density
# --------------------------------------------------------------------------
def mp_density(lam: np.ndarray, Q: float, sigma2: float = 1.0) -> np.ndarray:
    """
    Marchenko-Pastur probability density, for Q = T/N <= 1 (no point mass
    at zero in this regime).
    """
    lambda_min = sigma2 * (1 - np.sqrt(1/Q)) ** 2
    lambda_max = sigma2 * (1 + np.sqrt(1/Q)) ** 2

    density = np.zeros_like(lam)
    mask = (lam >= lambda_min) & (lam <= lambda_max)
    density[mask] = (
        np.sqrt((lambda_max - lam[mask]) * (lam[mask] - lambda_min))
        / (2 * np.pi * (1/Q) * sigma2 * lam[mask])
    )
    return density

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    n_species = 100
    scales = [7, 42]
    reference_census = 4
    n_spins = 1000
    standardize = True
    data_types = ('forest', 'senm')
    LABELS = {'forest': 'BCI', 'senm': 'SENM'}
    COLORS = {'forest': '#2ca02c', 'senm': '#1f77b4'}
    SENM_XMAX = {7: 6, 42: 19}   # scale -> right limit of the SENM panels
    rng = np.random.default_rng(0)

    for data_type in data_types:          # data_types = ['forest', 'senm'] (unique types, not one per panel)
        fig, axes = plt.subplots(1, len(scales), figsize=(6 * len(scales), 5))
        if len(scales) == 1:
            axes = [axes]

        for ax, scale in zip(axes, scales):
            print(f"Data type: {data_type}")
            print(f"Scale {scale}: building spatial series...")
            names = top_species(data_type, n_species, reference_census)
            df = raw2df(data_type, reference_census, names)
            spatial_series = compute_spatial_series(
                df, scale=scale, data_type=data_type, standardize=standardize
            )

            N, T = spatial_series.shape
            Q = T / N
            print(f"Scale {scale}: N={N}, T={T}, Q={Q:.4f}")
            
            mode = 'rows'
            print(f"Scale {scale}: running {n_spins} {mode} spins...")
            eigs_rows = spin_null_eigenvalues(spatial_series, n_spins=n_spins, rng=rng, mode=mode)
            lambda_rows_max = np.max(eigs_rows)

            """
            mode = 'columns'
            print(f"Scale {scale}: running {n_spins} {mode} spins...")
            eigs_columns = spin_null_eigenvalues(spatial_series, n_spins=n_spins, rng=rng, mode=mode)
            """

            mode = 'both'
            print(f"Scale {scale}: running {n_spins} {mode} spins...")
            eigs_both = spin_null_eigenvalues(spatial_series, n_spins=n_spins, rng=rng, mode=mode)

            # Real (unspun) eigenvalues -- only N=100 points, too few for a clean
            # histogram (density=True is still valid, just noisy with this few
            # samples per bin), so show as a rug plot instead of a histogram.
            C_real = np.corrcoef(spatial_series)
            real_eigs = np.linalg.eigvalsh(C_real)

            # Empirical histograms
            ax.hist(eigs_rows, bins=60, density=True, alpha=0.6,
                    color=COLORS[data_type], label=f'Row-spin null ({n_spins} spins, pooled), {LABELS[data_type]}')
            ax.axvline(lambda_rows_max, color=COLORS[data_type], ls='--', lw=1, alpha=0.7)
            #ax.hist(eigs_columns, bins=60, density=True, alpha=0.6,
            #        color='r', label=f'Column-spin null ({n_spins} spins, pooled), {data_type}')
            ax.hist(eigs_both, bins=60, density=True, alpha=0.5,
                    color='darkorange', label=f'Row+col spin null ({n_spins} spins, pooled), {LABELS[data_type]}')

            # Theoretical MP curve
            if Q >= 1:
                lam_grid = np.linspace(1e-6, eigs_both.max() * 1.05, 1000)
                mp_curve = mp_density(lam_grid, Q)
                ax.plot(lam_grid, mp_curve, color='crimson', lw=2, label='Theoretical MP density')

                lambda_min = (1 - np.sqrt(1/Q)) ** 2
                lambda_max = (1 + np.sqrt(1/Q)) ** 2
                ax.axvline(lambda_min, color='crimson', ls='--', lw=1, alpha=0.7)
                ax.axvline(lambda_max, color='crimson', ls='--', lw=1, alpha=0.7)
            else:
                ax.text(0.5, 0.5, 'Q<1 : outside MP range', transform=ax.transAxes,
                        ha='center', va='center', color='crimson')

            # Rug plot: one tick per real eigenvalue, near the bottom of the axes
            ymax = ax.get_ylim()[1]
            tick_height = 0.04 * ymax
            ax.vlines(real_eigs, 0, tick_height, color='black', lw=3, alpha=0.8,
                       label=f'Real eigenvalues (N={n_species})')

            ax.set_title(f'Scale = {scale}  (N={N}, T={T}, Q={Q:.3f})')
            ax.set_xlabel('Eigenvalue')
            ax.set_ylabel('Density')
            ax.legend(fontsize=9)
            if data_type == 'senm':
                ax.set_xlim(right=SENM_XMAX[scale])

        fig.suptitle('Spin-null empirical eigenvalue distribution vs. theoretical Marchenko-Pastur')
        fig.tight_layout()

        out_path = f'figures_supp/spin_null_vs_mp_{data_type}.svg'
        fig.savefig(out_path, dpi=150)
        print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
