"""
Supplementary information figures.

SI.1 - 100_species_rank.pdf
    Rank (by abundance) of the 100 most abundant species from census 4,
    tracked across all censuses.

SI.2 - mp_Q_validity.pdf
    Q = T / N_s plotted against scale, N_s fixed at 100, showing how
    sampling quality degrades (Q -> 1) as scale grows.

Run: python si_figures.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.utils import raw2df, top_species
from src.compute import get_dimensions


OUT_DIR = "figures_supp"
LABELS = {'forest': 'BCI', 'senm': 'SENM'}
COLORS = {'forest': '#2ca02c', 'senm': '#1f77b4'}

# ============================================================================
# SI.1 - species rank trajectories
# ============================================================================

def full_community_rank_over_censuses(censuses: list[int]) -> pd.DataFrame:
    """
    Rank every species in the whole community by abundance (1 = most
    abundant), per census.

    ASSUMPTION (please verify): raw2df('forest', census, None) -- or
    however your real API returns an unfiltered census -- returns every
    species, not just a pre-selected subset. Every other usage we've seen
    of raw2df always passed an explicit names list, so this call may need
    adjusting to match your actual signature.

    Returns
    -------
    pd.DataFrame, index=all species names (union across censuses),
    columns=censuses, values=rank
    """
    dfs = {census: raw2df('forest', census, None) for census in censuses}

    all_names = sorted(set().union(*(df['name'].unique() for df in dfs.values())))

    ranks = {}
    for census, df in dfs.items():
        abundance = df.groupby('name').size().reindex(all_names, fill_value=0)
        ranks[census] = abundance.rank(method='min', ascending=False)
    return pd.DataFrame(ranks)


def plot_species_rank(n_species: int, reference_census: int, censuses: list[int]):
    print("SI.1: computing species rank trajectories (whole community)...")
    top_names = top_species('forest', n_species, reference_census)
    rank_df = full_community_rank_over_censuses(censuses)

    # top_names are, by construction, the n_species most abundant at
    # reference_census -- so their full-community rank there is exactly
    # 1..n_species (ties aside), same scale the color gradient used before.
    ref_rank = rank_df.loc[top_names, reference_census]
    norm_rank = (ref_rank - 1) / (n_species - 1)

    fig, ax = plt.subplots(figsize=(8, 6))

    other_names = rank_df.index.difference(top_names)
    for name in other_names:
        ax.plot(censuses, rank_df.loc[name], color='lightgray', alpha=0.3, lw=0.6, zorder=1)

    cmap = plt.get_cmap('viridis')
    for name in top_names:
        ax.plot(censuses, rank_df.loc[name], color=cmap(norm_rank[name]),
                alpha=0.8, lw=1.2, zorder=2)

    ax.invert_yaxis()  # rank 1 at top
    ax.set_xlabel('Census')
    ax.set_ylabel('Rank in whole community (1 = most abundant)')
    ax.set_title(f'Rank of the {n_species} most abundant species '
                 f'(census {reference_census}) within the whole community')
    ax.set_xticks(censuses)

    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color=cmap(0.5), lw=1.5, label=f'Top {n_species} (census {reference_census})'),
        Line2D([0], [0], color='lightgray', lw=1.5, label='Rest of community'),
    ]
    ax.legend(handles=legend_elems, fontsize=9, loc='lower left')

    fig.tight_layout()
    out_path = f"{OUT_DIR}/100_species_rank.pdf"
    fig.savefig(out_path)
    print(f"Saved {out_path}")


# ============================================================================
# SI.2 - Q = T / N_s vs scale
# ============================================================================

def plot_mp_Q_validity(n_species: int, scales: list[int]):
    
    print("SI.2: computing Q = T/N_s vs scale...")
    
    fig, ax = plt.subplots(figsize=(7, 5))
    for data_type in ('forest', 'senm'):
        grid_width, grid_height = get_dimensions(data_type)

        T_values = []
        for scale in scales:
            n_bins_x = int(grid_width / scale)
            n_bins_y = int(grid_height / scale)
            T_values.append(n_bins_x * n_bins_y)
        Q_values = np.array(T_values) / n_species

        ax.semilogy(scales, Q_values, 'o-', label=LABELS[data_type], color=COLORS[data_type])

    ax.axhline(1.0, color='crimson', ls='--', lw=1, label='Q=1 (MP validity limit)')
    ax.set_xlabel('Scale (m)')
    ax.set_ylabel(r'$Q = T / N_s$')
    ax.set_title(f'Sampling quality vs. scale (N={n_species} species fixed)')
    ax.legend(fontsize=9)

    fig.tight_layout()
    out_path = f"{OUT_DIR}/mp_Q_validity.pdf"
    fig.savefig(out_path)
    print(f"Saved {out_path}")


def main():
    n_species = 100
    reference_census = 4
    censuses = list(range(1, 9))  # adjust to your actual census list / config['censuses']
    scales = list(range(5, 50))   # 5 to 33 inclusive; adjust step if not every integer

    plot_species_rank(n_species, reference_census, censuses)
    plot_mp_Q_validity(n_species, scales)


if __name__ == "__main__":
    main()
