import numpy as np
import os
import pandas as pd
from tqdm import tqdm

from src.config import load_config
from src.utils import get_dimensions

# === Load config ===
config = load_config()

def compute_spatial_series(df: pd.DataFrame, 
                           scale: int, data_type: str, standardize: bool = True):

    """Compute species abundance matrix.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input data with columns 'x', 'y', and a species identifier column.
    scale : int 
        length of box for construction of the series (spatial scale)
    data_type : str
        'forest' or 'senm'
    standardize : bool, optional
        If True, standardize each species’ abundance vector (z-score normalization).
    """

    grid_width, grid_height = get_dimensions(data_type)

    n_bins_x = int(grid_width/scale)
    n_bins_y = int(grid_height/scale)

    n_species = df['name'].nunique()

    df = df.copy()
    
    # Assign bins
    df['x_bin'] = (df['x'] / scale).astype(int).clip(0, n_bins_x - 1)
    df['y_bin'] = (df['y'] / scale).astype(int).clip(0, n_bins_y - 1)
    
    # Initialize spatial_series matrix
    spatial_series = np.zeros((n_species, n_bins_x * n_bins_y))
    
    for i, (_, group) in enumerate(df.groupby('name')):
        if i >= n_species:
            break
        bin_counts = group.groupby(['x_bin', 'y_bin']).size()
        for (x, y), count in bin_counts.items():
            spatial_series[i, x * n_bins_y + y] = count
    
    # Standardize each row individually if requested
    if standardize:
        row_means = spatial_series.mean(axis=1, keepdims=True)
        row_stds = spatial_series.std(axis=1, keepdims=True)
        # Avoid division by zero
        row_stds[row_stds == 0] = 1
        spatial_series = (spatial_series - row_means) / row_stds
    
    return spatial_series


def compute_forest_C(n_species: int, scale: int,
                     censuses: list[int] = config['censuses'],
                     reference_census: int = 4,
                     standardize: bool = True):
    """
    (docstring unchanged, add:)
    Returns
    -------
    ...
    spatial_series_list : list of np.ndarray
        Per-census spatial_series matrices, shape (n_species, n_bins) each.
    """
    from src.utils import raw2df, top_species

    names = top_species('forest', n_species, reference_census)
    matrices = []
    spatial_series_list = []

    for census in tqdm(censuses, desc='Computing forest C...'):
        df        = raw2df('forest', census, names)
        abundance = compute_spatial_series(df, scale=scale,
                                           data_type='forest',
                                           standardize=standardize)
        matrices.append(np.corrcoef(abundance))
        spatial_series_list.append(abundance)

    matrices = np.array(matrices)  # (n_censuses, n_species, n_species)
    return matrices.mean(axis=0), matrices.std(axis=0), names, spatial_series_list

def compute_senm_spectrum(n_species: int, scale: int,
                  n_censuses: int = 100, standardize: bool = True):
    """
    (docstring unchanged, add:)
    Returns
    -------
    ...
    spatial_series_list : list of np.ndarray
        Per-realization spatial_series matrices, shape (n_species, n_bins) each.
    """
    from src.utils import raw2df, top_species

    spectra = []
    spatial_series_list = []

    for census in tqdm(range(1, n_censuses + 1), desc='Computing senm spectra...'):
        names = top_species('senm', n_species, census)
        df    = raw2df('senm', census, names)

        abundance = compute_spatial_series(df, scale=scale,
                                           data_type='senm',
                                           standardize=standardize)

        C    = np.corrcoef(abundance)
        eigs = np.sort(np.linalg.eigvalsh(C))[::-1]  # descending
        spectra.append(eigs)
        spatial_series_list.append(abundance)

    spectra = np.array(spectra)  # (n_censuses, n_species)

    mean_spectrum = spectra.mean(axis=0)
    std_spectrum  = spectra.std(axis=0)

    return mean_spectrum, std_spectrum, spectra, spatial_series_list
 
def compute_forest_spectrum(n_species: int, scale: int,
                    censuses: list[int] = config['censuses'] , 
                    reference_census: int = 4, standardize: bool = True):
    """
    Compute mean eigenvalue spectrum for forest data across multiple censuses.

    Species are selected once from reference_census and tracked across all
    censuses.

    Parameters
    ----------
    n_species : int
        Number of most-abundant species to use, ranked from reference_census.
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
    mean_spectrum : np.ndarray, shape (n_species,)
        Mean eigenvalue spectrum, sorted descending, averaged across censuses.
    std_spectrum : np.ndarray, shape (n_species,)
        Standard deviation of each ranked eigenvalue across censuses.
    """
    from src.utils import raw2df, top_species

    names = top_species('forest', n_species, reference_census)

    spectra = []

    for census in tqdm(censuses, desc='Computing forest spectra...'):
        df        = raw2df('forest', census=census, names=names)
        abundance = compute_spatial_series(df, scale=scale,
                                           data_type='forest',
                                           standardize=standardize)

        C    = np.corrcoef(abundance)
        eigs = np.sort(np.linalg.eigvalsh(C))[::-1]  # descending
        spectra.append(eigs)

    spectra = np.array(spectra)  # (n_censuses, n_species)

    mean_spectrum = spectra.mean(axis=0)
    std_spectrum  = spectra.std(axis=0)

    return mean_spectrum, std_spectrum, spectra

def marchenko_pastur_bounds(data_type: str, n_species: int, scale: int):
    """Calculate the Marchenko-Pastur bounds."""

    grid_width, grid_height = get_dimensions(data_type)
    
    N = n_species
    T = (grid_height/scale) * (grid_width/scale)
    Q = T/N

    if Q < 1 :
        print('Out of MP range')
        lambda_min = 0
        lambda_max = 0
    else : 
        lambda_min = (1 - np.sqrt(1/Q))**2
        lambda_max = (1 + np.sqrt(1/Q))**2
     
    return lambda_max, lambda_min


def marchenko_pastur_filter(C: np.ndarray, scale: int, data_type: str, 
              remove_largest: bool=False, remove_small: bool=True) -> np.ndarray:
    """Uses Marchenko-Pastur law to remove noise.
    

        remove_largest (bool), optional
            If ``False``, all the eigenvectors associated to the
            significant eigenvalues will be used to reconstruct the
            de-noised empirical correlation matrix. If ``True``, the
            eigenvector associated to the largest eigenvalue (normally
            known as the ``market`` mode) is going to be excluded from
            the recontruction step. """ 
            

    n_species = np.size(C, 0) 

    grid_width, grid_height = get_dimensions(data_type)

    N = n_species
    T = (grid_height/scale) * (grid_width/scale)
    Q = T/N  
           
    if N > T:
        raise ValueError("T must be greater or equal than N.")

    w, v = np.linalg.eigh(C)  # Spectral decomposition of C
    

    w_max, w_min = marchenko_pastur_bounds(data_type, N, scale)
    
    if remove_small:
        selected = w > w_max
    else:
        selected = (w < w_min) | (w > w_max)
        

    if remove_largest:
        selected[-1] = False

    w_signal = w[selected]
    v_signal = v[:, selected]

    C_MP = v_signal.dot(np.diag(w_signal)).dot(v_signal.T)

    return C_MP


