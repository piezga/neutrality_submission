import numpy as np 
import networkx as nx 
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.linalg import expm

def L_detect_communities(corr_matrix, n_communities=None, tau=1e-3, Th=1e-4, return_linkage=False):
    """
    Detect communities from filtered correlation matrix using Laplacian diffusion clustering.
    
    Parameters
    ----------
    corr_matrix : array
        Correlation matrix
    n_communities : int, optional
        If provided, cuts dendrogram to get exactly this many communities.
        If None, uses Th threshold with distance criterion.
    tau : float
        Diffusion time parameter
    Th : float
        Distance threshold for cutting (only used if n_communities is None)
    return_linkage : bool
        If True, also return linkage matrix and cut height
        
    Returns
    -------
    reordered : array
        Correlation matrix reordered by community membership
    CM : array
        Community membership labels
    idx : array
        Indices for reordering
    linkage_matrix : array (optional)
        Normalized linkage matrix, returned if return_linkage=True
    cut_height : float (optional)
        Height at which dendrogram was cut, returned if return_linkage=True
    """
    corr_pos = np.copy(corr_matrix)
    corr_pos[corr_pos < 0] = 0  # only positive correlations
    
    # Build graph and Laplacian
    G = nx.from_numpy_array(np.abs(corr_pos))
    G.remove_edges_from(nx.selfloop_edges(G))
    L = nx.laplacian_matrix(G).todense()
    
    # Diffusion process
    num = expm(-tau * L)
    rho = num / np.trace(num)
    
    # Symmetric distance matrix
    Trho = np.copy(1.0 / rho)
    Trho = np.tril(Trho) + np.triu(Trho.T, 1)
    np.fill_diagonal(Trho, 0)
    
    # Hierarchical clustering
    dists = squareform(Trho)
    linkage_matrix = linkage(dists, "ward")
    
    # Normalize linkage distances by maximum
    linkage_matrix_norm = linkage_matrix.copy()
    linkage_matrix_norm[:, 2] = linkage_matrix[:, 2] / linkage_matrix[-1, 2]
    
    # Cut dendrogram
    if n_communities is not None:
        # Use maxclust criterion to get exactly n_communities
        CM = fcluster(linkage_matrix_norm, t=n_communities, criterion='maxclust')
        
        # Calculate the cut height for n_communities
        n_samples = linkage_matrix_norm.shape[0] + 1
        merge_index = n_samples - n_communities - 1
        cut_height = linkage_matrix_norm[merge_index, 2] if merge_index >= 0 else 0
    else:
        # Use distance threshold
        CM = fcluster(linkage_matrix_norm, t=Th, criterion="distance")
        cut_height = Th
    
    # Reorder matrix by community
    idx = np.argsort(CM)
    reordered = np.array([[corr_matrix[i][j] for j in idx] for i in idx])
    
    if return_linkage:
        return reordered, CM, idx, linkage_matrix_norm, cut_height
    else:
        return reordered, CM, idx
