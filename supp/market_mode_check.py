"""
Check whether the top eigenvector ("market mode" candidate) has uniform
sign (density mode) or a mixed-sign split (two opposing groups).

Panel 1: sorted loadings, colored by sign.
Panel 2: loadings vs. abundance rank (rules out "it's just density").

Run: python market_mode_check.py
"""

import numpy as np
import matplotlib.pyplot as plt

from src.compute import compute_forest_C


OUT_DIR = "figures_supp"

POS_COLOR = "#2E86AB"
NEG_COLOR = "#E07A5F"


def plot_market_mode_check(n_species: int, scale: int, censuses, reference_census: int):
    print(f"Computing forest C at scale {scale}m...")
    mean_C, std_C, names, spatial_series_list = compute_forest_C(
        n_species, scale, censuses, reference_census
    )

    w, v = np.linalg.eigh(mean_C)      # ascending
    top = v[:, -1]                     # eigenvector of the largest eigenvalue

    # names are assumed ordered by abundance at reference_census (rank 1 = index 0)
    rank = np.arange(1, n_species + 1)

    frac_pos = np.mean(top > 0)
    print(f"Top eigenvector: {frac_pos*100:.1f}% positive-signed loadings "
          f"({(1 - frac_pos)*100:.1f}% negative) -- "
          f"{'uniform (density-like)' if frac_pos in (0, 1) else 'mixed sign'}")

    colors = np.where(top > 0, POS_COLOR, NEG_COLOR)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Panel 1: sorted loadings, colored by sign ---
    ax = axes[0]
    order = np.argsort(top)
    ax.bar(np.arange(n_species), top[order], color=colors[order], width=0.9)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xlabel('Species (sorted by loading)')
    ax.set_ylabel('Loading')
    ax.set_title('Top eigenvector loadings')
    ax.spines[['top', 'right']].set_visible(False)

    # --- Panel 2: loadings vs. abundance rank ---
    ax = axes[1]
    ax.scatter(rank, top, c=colors, s=18, alpha=0.8, edgecolors='none')
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xlabel('Abundance rank')
    ax.set_ylabel('Loading')
    ax.set_title('Loading vs. abundance rank')
    ax.spines[['top', 'right']].set_visible(False)

    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=POS_COLOR, markersize=8, label='Positive'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=NEG_COLOR, markersize=8, label='Negative'),
    ]
    fig.legend(handles=legend_elems, loc='upper center', ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.04), fontsize=9)

    fig.suptitle(f'Top eigenvector sign structure at scale {scale}m', y=1.08)
    fig.tight_layout()

    out_path = f"{OUT_DIR}/market_mode_check.pdf"
    fig.savefig(out_path, bbox_inches='tight')
    print(f"Saved {out_path}")


def main():
    n_species = 100
    reference_census = 4
    censuses = list(range(1, 9))  # adjust to your actual census list
    scale = 40                    # adjust to your representative scale

    plot_market_mode_check(n_species, scale, censuses, reference_census)


if __name__ == "__main__":
    main()
