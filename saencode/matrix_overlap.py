# Import libraries
import numpy as np
from scipy.linalg import sqrtm


def check_mat_is_psd(mat: np.ndarray, atol: float=1e-8) -> bool:
    """
    True if M is symmetric and positive semi-definite.
    
    Parameters
    ----------
    mat : np.ndarray
        Symmetric normalized mutual information or covariance matrix (n, n).
    atol : float
        Absolute tolerance for comparing differences of two matrices.

    Returns
    -------
    bool
        True, if matrix is symmetric and positive semi-definite.
    """

    # Check that the matrix is 2-dimensional
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        return False
    # Check that matrix is symmetrical
    elif not np.allclose(mat, mat.T, atol=atol):
        return False
    # Check if matrix is positive semi-definite
    else:
        return bool(np.all(np.linalg.eigvalsh(mat) >= -atol))


def matrix_overlap(A: np.ndarray, B: np.ndarray) -> float:
    """ 
    Compute matrix overlap between two matrices.

    Implementation is based on the following reference:
    Berk Hess on "Convergence of sampling in protein simulations"
    Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 2002.
    
    Parameters
    ----------
    A : np.ndarray
        A symmetric matrix.
    B : np.ndarray
        A symmetric matrix.

    Returns
    -------
    overlap : float
        A value between [0, 1] describing similarity between two matrices.
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

def compute_matrix_overlap(nmi_blocks: list[np.ndarray]) -> np.ndarray:
    """ 
    Compute matrix overlap for a list of nMI blocks.
    
    Parameters
    ----------
    nmi_blocks : list
        A list of nMI blocks (matrices) for a given trajectory.
    
    Returns
    -------
    overlap : np.ndarray
        A matrix (num_of_blocks, num_of_blocks) describing pairwise overlap between nMI blocks for a given trajectory.

    """

    overlap = np.zeros((len(nmi_blocks), len(nmi_blocks)))
    for i in range(len(nmi_blocks)):
        for j in range(len(nmi_blocks)):
            overlap[i, j] = matrix_overlap(nmi_blocks[i], nmi_blocks[j])

    return overlap
