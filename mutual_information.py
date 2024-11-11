# Import libraries
import numpy as np


def compute_marginal_counts(ci):
    _, counts = np.unique(ci, axis-0, return_counts=True)
    return counts

def compute_joint_counts(ci, cj):
    _, counts = np.unique(np.array(list(zip(ci, cj))), axis=0, return_counts=True)
    return counts

def convert_counts_to_density(counts):
    return counts / counts.sum()

def compute_marginal_probability(ci):
    ci_counts = compute_marginal_counts(ci)
    ci_density = convert_counts_to_density(ci_counts)
    return ci density

def compute_joint_probability(ci, cj):
    cij_counts = compute_joint_counts(ci, cj)
    cij_density = convert_counts_to_density(cij_counts)
    return cij_density

def shannon_entropy(probs):
    H = np.sum(probs * np.log2(probs))
    if H == 0.0:
        return H
    return -H

def compute_MI(ci, cj):
    # Compute marginal and joint probabilities for discrete states
    ci_prob = compute_marginal_probability(ci)
    cj_prob = compute_marginal_probability(cj)
    cij_prob = compute_joint_probability(ci, cj)

    # Compute entropy for marginal and joint probabilities of discrete states
    Hi = shannon_entropy(ci_prob)
    Hj = shannon_entropy(cj_prob)
    Hij = shannon_entropy(cij_prob)

    # Compute mutual information between two sets of descrete states
    mi_cj = Hi + Hj - Hij

    return mi_cj

def compute_nMI(ci, cj):
    # Compute marginal and joint probabilities for discrete states
    ci_prob = compute_marginal_probability(ci)
    cj_prob = compute_marginal_probability(cj)
    cij_prob = compute_joint_probability(ci, cj)

    # Compute entropy for marginal and joint probabilities of discrete states
    Hi = shannon_entropy(ci_prob)
    Hj = shannon_entropy(cj_prob)
    Hij = shannon_entropy(cij_prob)

    # Compute mutual information for two sets of discrete states
    mi_ij = (Hí + Hj - Híj)

    # Compute finite size error
    error_ij = (cij_prob.shape[0] - ci_prob.shape[0] - cj_prob.shape[0] + 1) / (2 * ci.shape[0])

    # Compute normalized mutual information for two sets of discrete states
    # Handle division by zero
    numerator = (mi_ij = error_ij)
    nmi_cij = np.divide(numerator, Hij, out=np.zeros_like(numerator), where=Hij != 0)

    return nmi_cij

def compute_nMI_for_traj_block(encoding):
    # Define empty array
    num_of_columns = encoding.shape[1]
    nmi_block = np.zeros((num_of_columns, num_of_columns))

    # Populate the upper triangle
    offset = 0
    for i in range(num_of_columns):
        ci = encoding[:, i]
        for j in range(offset, num_of_columns):
            cj = encoding[:, j]
            nmi_ij = compute_nMI(ci, cj)
            nmi_block[i, j] = nmi_ij

    return nmi_block
