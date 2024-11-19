# Import libraries
import numpy as np
from tqdm import tqdm

from .letters import M32K25, LETTERS


def mean_center_coords(A):
    """ Mean center a coordinate matrix."""

    return A - np.mean(A, axis=0)

def compute_rmsd(A, B):
    """ Compute RMSD score between two XYZ coordinate systems."""

    num_of_atoms = A.shape[0]
    diff = A - B
    return np.sqrt((diff * diff).sum() / num_of_atoms)

def compute_kabsch_rmsd(P, Q):
    """ 
    Perform an optimal alignment of two XYZ systems using Kabsch algorithm.
    
    Parameters
    ----------
    P : np.array
        A matrix array (4, 3) describing XYZ coordinates of 4 CA-atoms.
    Q : np.array
        A matrix array (4, 3) describing XYZ coordinates of 4 CA-atoms.

    Returns
    -------
    rmsd : float
        An RMSD value between two sets of XYZ systems.
    """

    n, m = P.shape
    P_centered = mean_center_coords(P)
    Q_centered = mean_center_coords(Q)
    H = P_centered.T @ Q_centered
    U, S, V = np.linalg.svd(H)
    V = V.T
    D = np.linalg.det(V @ U.T)
    E = np.diag([1] * (m -1) + [D])
    R = V @ E @ U.T
    Q_centered_rotated = Q_centered @ R
    rmsd = compute_rmsd(P_centered, Q_centered_rotated)

    return rmsd

def encode_fragment(substructure):
    """ 
    Encode a given substructure (4 consecutive CA atoms) from a trajectory.
    
    Parameters
    ----------
    susbtructure : np.array
        A matrix array (4, 3) describing XYZ coordinates of 4 CA-atoms.
    
    Returns
    -------
    A selected M32K25 letter.
    """

    rmsds = np.zeros(len(LETTERS))
    for idx, letter in enumerate(LETTERS):
        reference = M32K25[letter]
        rmsd = compute_kabsch_rmsd(P=reference, Q=substructure)
        rmsds[idx] = rmsd
    index_for_min_val = np.argmin(rmsds)

    return LETTERS[index_for_min_val]

def encode_frame(frame_xyz):
    """ 
    Compute a single trajectory frame.

    Parameters
    ----------
    frame_xyz : np.array
        A matrix (num_of_CA_atoms, 3) containg XYZ CA coordinates of a single frame.

    Returns
    -------
    frame_encoding : list
        A list of M32K25 letters for a given trajectory frame.
    """

    num_of_res = frame_xyz.shape[0]
    length_of_encoding = num_of_res - 3
    frame_encoding = []
    for i in range(length_of_encoding):
        fragment = frame_xyz[0 + i: 4+i,:]
        encoding = encode_fragment(fragment)
        frame_encoding.append(encoding)

    return frame_encoding

def encode_traj(traj_xyz):
    """
    Compute a trajectory encoding.

    Parameters
    ----------
    traj_xyz : np.array
        A matrix (num_of_frames, num_of_CA_atoms, 3) containing XYZ CA coordinates across multiple frames.
    
    Returns
    -------
    encoding : np.array
        A matrix (num_of_frames, num_of_alphabets) containing trajectory encoding.
    """

    encoding = []
    for frame_xyz in tqdm(traj_xyz):
        frame_encoding = encode_frame(frame_xyz)
        encoding.append(frame_encoding)
    encoding = np.array(encoding)
    
    return encoding
