# Import libraries
import numpy as np
from scipy.linalg import sqrtm


def check_mat_is_psd(mat, tol=1e-8):
    """
    True if M is symmetric and positive semi-definite.
    
    Parameters
    ----------
    mat : np.array
        Symmetric normalized mutual information or covariance matrix (n, n).

    Returns
    -------
    bool
        True, if matrix is symmetrix and positive semi-definite.
    """

    # Check that the matrix is 2-dimensional
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        return False
    # Check that matrix is symmetrical
    elif not np.allclose(mat, mat.T, atol=tol):
        return False
    # Check if matrix is positive semi-definite
    else:
        return bool(np.all(np.linalg.eigvalsh(mat) >= -tol))


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

    # Check if matrices are positive semi-definite
    if not (check_mat_is_psd(A) and check_mat_is_psd(B)):
        raise ValueError("matrix_overlap requires symmetric PSD matrices")

    # Compute overlap between two matrices
    diff = sqrtm(A) - sqrtm(B)
    numerator = np.sqrt(np.trace(diff @ diff))
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
        A matrix (num_of_blocks, num_of_blocks) describing pairwise overlap between nMI blocks for a given trajectory.

    """

    overlap = np.zeros((len(nmi_blocks), len(nmi_blocks)))
    for i in range(len(nmi_blocks)):
        for j in range(len(nmi_blocks)):
            overlap[i, j] = matrix_overlap(nmi_blocks[i], nmi_blocks[j])

    return overlap
