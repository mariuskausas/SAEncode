# Import libraries
import numpy as np
from tqdm import tqdm


def compute_marginal_counts(ci: np.ndarray) -> np.ndarray:
    """
    Count letters for a given column of encoding.

    Parameters
    ----------
    ci : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.

    Returns
    -------
    counts : np.ndarray
        A matrix (num_of_letters,) of respective counts for each letter.
    """

    _, counts = np.unique(ci, axis=0, return_counts=True)
    return counts


def compute_joint_counts(ci: np.ndarray, cj: np.ndarray) -> np.ndarray:
    """
    Count joint letters for two given columns of encoding.

    Parameters
    ----------
    ci : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.
    cj : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.

    Returns
    -------
    counts : np.ndarray
        A matrix (num_of_letter_pairs,) of respective counts for each letter pair.
    """

    _, counts = np.unique(np.array(list(zip(ci, cj))), axis=0, return_counts=True)
    return counts


def convert_counts_to_density(counts: np.ndarray) -> np.ndarray:
    """ Convert count data into a probability distribution."""
    return counts / counts.sum()


def compute_marginal_probability(ci: np.ndarray) -> np.ndarray:
    """
    Compute marginal letter probability for a single column in encoding.

    Parameters
    ----------
    ci : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.

    Returns
    -------
    ci_density : np.ndarray
        A matrix (num_of_letters,) of respective letter probabilities for a given column.
    """

    ci_counts = compute_marginal_counts(ci)
    ci_density = convert_counts_to_density(ci_counts)
    return ci_density


def compute_joint_probability(ci: np.ndarray, cj: np.ndarray) -> np.ndarray:
    """
    Compute joint letter probability for two columns in encoding.

    Parameters
    ----------
    ci : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.
    cj : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.

    Returns
    -------
    cij_density : np.ndarray
        A matrix (num_of_letter_pairs,) of respective letter-pair probabilities for two columns.
    """
    cij_counts = compute_joint_counts(ci, cj)
    cij_density = convert_counts_to_density(cij_counts)
    return cij_density


def shannon_entropy(probs: np.ndarray) -> float:
    """ Compute Shannon Entropy for a given probability distribution."""

    H = np.sum(probs * np.log2(probs))
    if H == 0.0:
        return H
    return -H


def compute_nMI(ci: np.ndarray, cj: np.ndarray) -> float:
    """
    Compute normalized mutual information (nMI) between two sets of encoding columns.

    Parameters
    ----------
    ci : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.
    cj : np.ndarray
        A matrix (num_of_trajectory_frames,) containing a column of encoding letters.
    
    Returns
    -------
    nmi_cij : float
        A nMI for two sets of encoding columns.
    """


    # Compute marginal and joint probabilities for discrete states
    ci_prob = compute_marginal_probability(ci)
    cj_prob = compute_marginal_probability(cj)
    cij_prob = compute_joint_probability(ci, cj)

    # Compute entropy for marginal and joint probabilities of discrete states
    Hi = shannon_entropy(ci_prob)
    Hj = shannon_entropy(cj_prob)
    Hij = shannon_entropy(cij_prob)

    # Compute mutual information for two sets of discrete states
    mi_ij = (Hi + Hj - Hij)

    # Compute finite size error
    error_ij = (cij_prob.shape[0] - ci_prob.shape[0] - cj_prob.shape[0] + 1) / (2 * ci.shape[0])

    # Compute normalized mutual information for two sets of discrete states
    # Handle division by zero
    numerator = (mi_ij - error_ij)
    nmi_cij = np.divide(numerator, Hij, out=np.zeros_like(numerator), where=Hij != 0)

    return float(nmi_cij)


def compute_nMI_for_single_block(encoding: np.ndarray) -> np.ndarray:
    """ 
    Compute a normalized mutual information (nMI) for a single encoding block.

    Parameters
    ----------
    encoding : np.ndarray
        A matrix (num_of_frames, num_of_alphabet_residues) containing trajectory encoding.

    Returns
    -------
    nmi_block : np.ndarray
        A matrix (num_of_alphabet_residues, num_of_alphabet_residues) containing nMI between each encoding position.

    """

    # Define empty array
    num_of_columns = encoding.shape[1]
    nmi_block = np.zeros((num_of_columns, num_of_columns))

    # Compute nMI
    for i in range(num_of_columns):
        ci = encoding[:, i]
        for j in range(num_of_columns):
            cj = encoding[:, j]
            nmi_ij = compute_nMI(ci, cj)
            nmi_block[i, j] = nmi_ij
            
    return nmi_block


def compute_nMI_for_traj_blocks(encoding: np.ndarray, number_of_blocks: int = 5) -> list[np.ndarray]:
    """ 
    Compute a normalized mutual information (nMI) for multiple encoding blocks.

    Parameters
    ----------
    encoding : np.ndarray
        A matrix (num_of_frames, num_of_alphabet_residues) containing trajectory encoding.
    number_of_blocks : int
        Number of blocks to cut trajectory into.

    Returns
    -------
    traj_nmi_blocks : list[np.ndarray]
        A list of nMI matrices, one per block, each of shape
        (num_of_alphabet_residues, num_of_alphabet_residues).

    """

    # Calculate block size
    traj_nmi_blocks = []
    block_size = int(encoding.shape[0] / number_of_blocks)

    # Enumerate each block and compute nMI for a single block
    for block_id in tqdm(range(number_of_blocks)):
        encoding_per_block = encoding[block_size * block_id: block_size * (block_id + 1)]
        traj_nmi_per_block = compute_nMI_for_single_block(encoding=encoding_per_block)
        traj_nmi_blocks.append(traj_nmi_per_block)

    return traj_nmi_blocks
