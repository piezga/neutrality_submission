"""
Produces one plot:
  1. Scatter plot of raw correlations vs nutrient-regressed correlations after MP
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

from src.config import load_config
from src.compute import (compute_spatial_series, compute_forest_C, marchenko_pastur_filter)
from src.utils import raw2df
from tqdm import tqdm

# ============================================================
# Config
# ============================================================
config      = load_config()
n_species   = config['n_species']
scale       = config['scale']

output_dir = 'figures/panel_2/'
os.makedirs(output_dir, exist_ok=True)

# ============================================================
# Functions 
# ============================================================

def load_nutrient_data(scale: int = 20, as_dataframe: bool = False, standardize: bool = True):
    """
    Load and bin BCI soil nutrient data to a given spatial scale.

    Parameters
    ----------
    scale : int
        Spatial binning scale in metres (default: 20).
    as_dataframe : bool
        If True, return a (n_sites, n_nutrients) DataFrame with column names.
        If False, return a (n_nutrients, n_sites) numpy array.
    standardize : bool
        If True, z-score each nutrient across sites.

    Returns
    -------
    np.ndarray of shape (n_nutrients, n_sites)  if as_dataframe=False
    pd.DataFrame of shape (n_sites, n_nutrients) if as_dataframe=True
    """
    from src.utils import get_dimensions

    nutrient_columns = ['Al', 'B', 'Ca', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'P', 'Zn', 'N', 'N(min)', 'pH']
    nutrient_file    = 'data/soil/barro_soil_data.xls'
    grid_width, grid_height = get_dimensions('forest') 

    n_bins_x = int(grid_width  / scale)
    n_bins_y = int(grid_height / scale)
    n_sites  = n_bins_x * n_bins_y

    df = pd.read_excel(nutrient_file)
    df['x_bin'] = (df['x'] / scale).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / scale).astype(int).clip(0, n_bins_y - 1)

    binned = df.groupby(['x_bin', 'y_bin'])[nutrient_columns].mean()

    # Build (n_nutrients, n_sites) array, NaN where no soil data
    nutrient_data = np.full((len(nutrient_columns), n_sites), np.nan)
    for x_bin in range(n_bins_x):
        for y_bin in range(n_bins_y):
            idx = x_bin * n_bins_y + y_bin
            if (x_bin, y_bin) in binned.index:
                nutrient_data[:, idx] = binned.loc[(x_bin, y_bin)].values

    # Fill missing bins with per-nutrient mean
    for i in range(len(nutrient_columns)):
        col_mean = np.nanmean(nutrient_data[i])
        nutrient_data[i, np.isnan(nutrient_data[i])] = col_mean

    if standardize:
        means = nutrient_data.mean(axis=1, keepdims=True)
        stds  = nutrient_data.std(axis=1,  keepdims=True)
        stds[stds == 0] = 1
        nutrient_data = (nutrient_data - means) / stds

    if as_dataframe:
        return pd.DataFrame(nutrient_data.T, columns=nutrient_columns)

    return nutrient_data

def mean_forest_series(n_species: int, scale: int,
                     censuses: list[int] = config['censuses'],
                     reference_census: int = 4,
                     standardize: bool = True):
    """
    Compute mean spatial series for forest data across multiple censuses.

    Species are selected once from reference_census and tracked across all
    censuses.

    Parameters
    ----------
    n_species : int
        Number of most-abundant species, ranked from reference_census.
    scale : int
        Spatial binning scale in metres.
    censuses : list[int]
        Census indices to average over.
    reference_census : int
        Census from which to rank and select the top n_species.
    standardize : bool
        Whether to z-score each species' abundance vector before correlating.

    Returns
    -------
    mean_C : np.ndarray, shape (n_species, n_species)
        Mean correlation matrix averaged across censuses.
    std_C : np.ndarray, shape (n_species, n_species)
        Element-wise standard deviation across censuses.
    """
    from src.utils import raw2df, top_species

    names = top_species('forest', n_species, reference_census)
    matrices = []

    for census in tqdm(censuses, desc='Computing forest C...'):
        df        = raw2df('forest', census, names)
        abundance = compute_spatial_series(df, scale=scale,
                                           data_type='forest',
                                           standardize=standardize)
        matrices.append(abundance)

    matrices = np.array(matrices)  # (n_censuses, n_species, n_species)
    return matrices.mean(axis=0), matrices.std(axis=0)

def compute_regression(n_species: int, scale: int = 20):
    """
    Remove nutrient effects from species abundance, returning cleaned time series.
    Fits a linear + quadratic (squared) model for each nutrient.
    
    Parameters
    ----------
    mean_abundance : array, shape (n_species, n_sites)
        Species abundance data (rows = species, cols = spatial bins)
    nutrient_data : array, shape (n_nutrients, n_sites)
        Nutrient abundance data (rows = nutrients, cols = spatial bins)
        
    Returns
    -------
    regressed_corr : array, shape (n_species, n_species)
        Correlation matrix of regressed abundances
    regressed_abundance : array, shape (n_species, n_sites)
        Species abundance with nutrient effects removed (residuals)
    """
    from sklearn.linear_model import LinearRegression
    
    mean_abundance, _ = mean_forest_series(n_species, scale)
    n_species = mean_abundance.shape[0]
    nutrient_data = load_nutrient_data()

    X_species  = mean_abundance.T   # (n_sites, n_species)
    X_nutrients = nutrient_data.T      # (n_sites, n_nutrients)

    # Build design matrix: [N, N²] for each nutrient
    X_design = np.hstack([X_nutrients, X_nutrients ** 2])  # (n_sites, 2*n_nutrients)


    residuals = np.zeros_like(X_species)

    print("\nComputing regressed abundances (removing nutrient effects, linear + quadratic)...")
    print(f"  Regressing out {X_nutrients.shape[1]} nutrients "
          f"({X_design.shape[1]} terms) from {n_species} species...")

    for i in range(n_species):
        if (i + 1) % 20 == 0:
            print(f"    Processed {i+1}/{n_species} species...")

        reg = LinearRegression()
        reg.fit(X_design, X_species[:, i])
        residuals[:, i] = X_species[:, i] - reg.predict(X_design)

    regressed_abundance = residuals.T
    regressed_corr      = np.corrcoef(residuals.T)

    print(f"  ✓ regressed abundance shape: {regressed_abundance.shape}")
    print(f"  ✓ regressed correlation matrix computed")

    return regressed_corr, regressed_abundance


def main():

    # ============================================================
    # Calculate data
    # ============================================================
    C, _, names, _ = compute_forest_C(n_species, scale)
    C_MP = marchenko_pastur_filter(C, scale, data_type='forest')
    C_reg, _ = compute_regression(n_species, scale)
    C_reg_MP = marchenko_pastur_filter(C_reg, scale, data_type='forest')

# Scatter plot
    fig, ax = plt.subplots(figsize=(3.35, 3.35), dpi=300)

    triu_indices = np.triu_indices_from(C_MP, k=1)
    raw_upper = C_MP[triu_indices]
    regressed_upper = C_reg_MP[triu_indices]
    valid_mask = ~(np.isnan(raw_upper) | np.isnan(regressed_upper))
    raw_upper = raw_upper[valid_mask]
    regressed_upper = regressed_upper[valid_mask]

    ax.scatter(raw_upper, regressed_upper, alpha=0.8, s=1, c='#000000', 
               edgecolors='none', rasterized=True)



# existing scatter
    ax.scatter(raw_upper, regressed_upper, alpha=0.8, s=1, c='#000000', 
               edgecolors='none', rasterized=True)

    pairs = [(42,74), (15, 41), (68, 30)]

    for i, j in pairs:
        ax.scatter(C_MP[i, j], C_reg_MP[i, j], 
                   alpha=1.0, s=5, c='red', edgecolors='black', linewidth=0.5)
                   

# limits
    lim = max(abs(raw_upper.min()), abs(raw_upper.max()), 
              abs(regressed_upper.min()), abs(regressed_upper.max()))

    x_vals = np.linspace(-lim, lim, 500)

# soil mediation zone between y=0 and y=x
    ax.fill_between(x_vals, 0, x_vals, 
                    color='#4C72B0', alpha=0.15, label='Prevalence of soil mediation')

# ensure symmetric limits
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    ax.legend()
    ax.set_xlabel('MP-filtered correlation', fontsize=10)
    ax.set_ylabel('MP-filtered correlation after regression', fontsize=10)
    ax.tick_params(labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_aspect('equal')
    plt.savefig(f'{output_dir}scatter.svg')

# Species raster plots
# Params 

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
    plt.rcParams['font.size'] = 7
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['xtick.major.size'] = 2
    plt.rcParams['ytick.major.size'] = 2
# Background

    bci_img = plt.imread('data/soil/level.png')
# plots

    species_colors = [
        ['#E69F00', '#7A1F5C'],
        ['#1f77b4', '#8c564b'],
        ['#d62728', '#2ca02c'],
    ]

    panel_labels = ['b', 'c', 'd']
    panel_names  = ['panel_b_species', 'panel_c_species', 'panel_d_species']

    for i, (pair, colors, label, name) in enumerate(
            zip(pairs, species_colors, panel_labels, panel_names)):

        fig, ax = plt.subplots(figsize=(3.5, 2.0), dpi=300)
        ax.imshow(
            bci_img[::-1], origin='lower', alpha=0.9,
            extent=[0, 1000, 0, 500], aspect='equal'
        )
        for idx, sp_idx in enumerate(pair):
            sp_name = names[sp_idx]
            sp_data = raw2df('forest', 4, sp_name)
            #sp_data = df[df['name'] == sp_name.strip()]
            ax.scatter(
                sp_data['x'], sp_data['y'],
                s=0.2, color=colors[idx], alpha=0.9,
                rasterized=True, label=f'{sp_name[:15]}'
            )
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 500)
        ax.set_xlabel('X (m)', fontsize=7)
        ax.set_ylabel('Y (m)', fontsize=7)
        ax.legend(fontsize=5, frameon=False, loc='upper right',
                  markerscale=4, handletextpad=0.3)
        fig.tight_layout()
        fig.savefig(f'{output_dir}{name}.svg', format='svg', bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ {name}.svg")

if __name__ == "__main__":
    main()
