# Import libraries
import numpy as np

from saencode.letters import M32K25, LETTERS
from saencode.encode import encode_fragment, encode_frame, encode_traj
from saencode.mutual_information import compute_nMI
from saencode.matrix_overlap import matrix_overlap, compute_matrix_overlap


def test_reference_fragment_encodes_to_itself():
    """ Test that each fragment encodes itself."""

    for letter in LETTERS:
        assert encode_fragment(M32K25[letter]) == letter


def test_encoding_is_invariant_to_rotation_and_translation():
    """ Test that fragment encoding is invariant to rotation and translation."""

    # Define a rotation (along X-axis)
    theta = np.radians(90)
    c, s = np.cos(theta), np.sin(theta)
    rotation = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    # Define a translation 
    translation = np.array([2.0, 4.5, -3.0])
    
    for letter in LETTERS:
        moved = M32K25[letter] @ rotation.T + translation
        assert encode_fragment(moved) == letter


def test_encode_frame_shape():
    """ Test that encoding process returns the right frame encoding size len(num_of_alphabets)."""

    num_of_res = 10
    rng = np.random.default_rng(0)
    frame = rng.standard_normal((num_of_res, 3))
    print(frame)
    assert len(encode_frame(frame)) == num_of_res - 3


def test_encode_traj_shape():
    """ Test that encoding process returns the right trajectory encoding shape (num_of_frames, num_of_alphabets)."""

    num_of_res = 10
    num_of_frames = 5
    rng = np.random.default_rng(0)
    traj = rng.standard_normal((num_of_frames, num_of_res, 3))
    encoding = encode_traj(traj)
    assert encoding.shape == (num_of_frames, num_of_res - 3)


def test_nMI_is_order_independent():
    """ Test that normalized mutual is not dependend on the order on input of two different inputs."""

    rng = np.random.default_rng(0)
    ci = rng.choice(list("ABC"), size=500)
    cj = rng.choice(list("ABC"), size=500)
    assert np.isclose(compute_nMI(ci, cj), compute_nMI(cj, ci))


def test_overlap_for_identical_matrices():
    """
    Test that overlap of two identical matrices is equal to 1.
    """

    A = np.array([[2.0, 0.5], [0.5, 1.0]])
    assert np.isclose(matrix_overlap(A, A), 1.0)


def test_overlap_is_order_independent():
    """ Test that matrix overlap is not dependend on the order on input of two different matrices."""

    A = np.array([[2.0, 0.5], [0.5, 1.0]])
    B = np.eye(2)
    assert np.isclose(matrix_overlap(A, B), matrix_overlap(B, A))


def test_overlap_shape_and_diagonal():
    """ Test that overlap is a square matrix with a diagonal of 1. """

    A = np.array([[2.0, 0.5], [0.5, 1.0]])
    B = np.eye(2)
    overlap = compute_matrix_overlap([A, B])

    assert overlap.shape == (2, 2)
    assert np.allclose(np.diag(overlap), 1.0)
