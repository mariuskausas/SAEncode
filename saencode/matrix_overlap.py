# Import libraries
import numpy as np
from scipy.linalg import sqrtm


def matrix_overlap(A, B):
    """ 
    Compute matrix overlap between two matrices.

    Implementation as follows the following reference:
    Berk Hess on "Convergence of sampling in protein simulations"
    Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 2002.
    
    Parameters
    ----------
    A : np.array
        A symmetric matrix.
    B : np.array
        A symmetric matrix.

    Returns
    -------
    overlap : float
        A value between [0, 1] describing similarity between two nmatrices.
    """

    diff = sqrtm(A) - sqrtm(B)
    numerator = np.sqrt(np.trace(diff * diff))
    denominator = np.sqrt(np.trace(A) + np.trace(B))
    overlap = 1 - (numerator / denominator)

    return overlap

def compute_matrix_overlap(nmi_blocks):
    """ 
    Compute matrix overlap for a list of nMI blocks.
    
    Paramaters
    ----------
    nmi_blocks : list
        A list of nMI blocks (matrices) for a given trajectory.
    
    Returns
    -------
    overlap : np.array
        A matrix describing pairwise overlap between nMI blocks for a given trajectory.

    """

    overlap = np.zeros((len(nmi_blocks), len(nmi_blocks)))
    for i in range(len(nmi_blocks)):
        for j in range(len(nmi_blocks)):
            overlap[i, j] = matrix_overlap(nmi_blocks[i], nmi_blocks[j])

    return overlap
