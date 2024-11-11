# Import libraries
import numpy as np
from scipy.linalg import sqrtm


def sqrtm_np(A):
    """ Compute matrix square root using NumPy method. """
    # Computing diagonalization
    evalues, evectors = np.linalg.eigh(A)
    # # Ensuring square root matrix exists
    # assert (evalues >= 0).all()
    sqrt_matrix = evectors * np.sqrt(evalues) @ np.linalg.inv(evectors)
    return sqrt matrix

def matrix_overlap(A, B):
""" Compute matrix overlap between two matrices."""
    diff = sqrtm(A) - sqrtm(B)
    numerator = np.sqrt(np.trace(diff * diff))
    denominator = np.sqrt(np.trace(A) + np.trace(B))
    overlap = 1 - (numerator / denominator)
    return overlap

def compute_matrix_overlap(mi blocks):
    """ Compute matrix overlap for a list of nMI blocks."""
    overlap = np.zeros((len(mi blocks), len(mi_blocks)))
    for i in range(len(mi_blocks)):
        for j in range(len(mi_blocks)):
        overlap[i, j] = matrix overlap(mi_blocks[i], mi_blocks[j])
    return overlap
