import numpy as np 
import pandas as pd 
from src.config import load_config

config = load_config()

def raw2df(data_type: str, census: int, names: list = None) -> pd.DataFrame:
    """
    Load forest or SENM spatial data into a unified DataFrame.

    Parameters
    ----------
    data_type : str
        'forest' or 'senm'
    census : int
        Census number (forest) or realization index (senm), 1-based.
    names: list
        list of species to keep

    Returns
    -------
    pd.DataFrame
        Columns: ['x', 'y', 'name']
    """
    if data_type == 'forest':
        path = f"data/raw/bci.tree{census}.csv"
        try:
            df = pd.read_csv(path)
            df = df[['x', 'y', 'name']]
        except FileNotFoundError:
            print(f"Warning: forest file not found at '{path}'. Returning empty DataFrame.")
            return pd.DataFrame(columns=['x', 'y', 'name'])
        except Exception as e:
            print(f"Warning: error loading forest file '{path}': {e}. Returning empty DataFrame.")
            return pd.DataFrame(columns=['x', 'y', 'name'])

    elif data_type == 'senm':
        senm_cfg  = config['senm']
        kernel    = senm_cfg['kernel']
        nu        = f"{float(senm_cfg['nu']):e}"   # as v2d writes it in the file name
        nx        = senm_cfg['grid_width']
        ny        = senm_cfg['grid_height']

        path = (f"data/senm/"
                f"AConf_NX{nx}_NY{ny}_NU{nu}_Square_Kernsize{kernel}_rel.{census}")
        try:
            raw = np.loadtxt(path)
            df  = pd.DataFrame(raw, columns=['x', 'y', 'species_id'])
            df['name'] = 'species_' + df['species_id'].astype(int).astype(str)
            df = df[['x', 'y', 'name']]
        except FileNotFoundError:
            print(f"Warning: SENM file not found at '{path}'. Returning empty DataFrame.")
            return pd.DataFrame(columns=['x', 'y', 'name'])
        except Exception as e:
            print(f"Warning: error loading SENM file '{path}': {e}. Returning empty DataFrame.")
            return pd.DataFrame(columns=['x', 'y', 'name'])

    else:
        raise ValueError(f"Unknown data_type '{data_type}'. Must be 'forest' or 'senm'.")
    
    df = df.dropna()
    # Only keep select species

    if names is not None:
        if isinstance(names, str):
            names = [names]
        df = df[df['name'].isin(names)]

    return df.reset_index(drop=True)

def top_species(data_type: str, n_species: int, census: int) -> list[str]:
    """Return names of the n most abundant species in a given census."""
    df = raw2df(data_type, census)
    counts = df.groupby('name').size()

    if data_type not in ('forest', 'senm'):
        raise ValueError(f"Unknown data type. Must be either 'forest' or 'senm'")

    return counts.nlargest(n_species).index.tolist()

def get_dimensions(data_type: str):
    """Returns width and height of selected dataset"""

    if data_type == 'forest':
        forest_cfg = config['forest']
        grid_width = forest_cfg['grid_width']
        grid_height = forest_cfg['grid_height']

    elif data_type == 'senm':
        senm_cfg = config['senm']
        grid_width = senm_cfg['grid_width']
        grid_height = senm_cfg['grid_height']
        
    else: 
        raise ValueError(f"Unknown data type. Must be either 'forest' or 'senm'.")

    return grid_width, grid_height


def load_spectra(data_type: str, scale: int):
    """Load precomputed spectra for a given data type and scale."""
    base = f"data/processed/scale_{scale}/{data_type}"
    return (
        np.load(f"{base}_spectra.npy")
    )

def load_C(data_type: str, scale: int):
    """Load precomputed mean and std correlation matrix for a given data type and scale."""
    base = f"data/processed/scale_{scale}/{data_type}"
    return (
        np.load(f"{base}_correlation_mean.npy"),
        np.load(f"{base}_correlation_std.npy"),
    )

def ell_N(scales, bci_counts, senm_counts, tol=2):
    """Smallest scale from which the BCI mean stays within tol SD of the neutral ensemble at all larger scales.
    bci_counts: (n_scales, 8); senm_counts: (n_scales, 100)."""
    ok = np.abs(bci_counts.mean(1) - senm_counts.mean(1)) <= tol * senm_counts.std(1, ddof=1)
    for j in range(len(scales)):
        if ok[j:].all():
            return int(scales[j])
