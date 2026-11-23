## Structural encoding and mutual information analysis of molecular dynamics trajectories.

### Reference

Python implementation of structural encoding using M32K25 and mutual information analysis of molecular dynamics trajectories based on the following references:

- [Detection of allosteric signal transmission by information-theoretic analysis of protein dynamics](https://faseb.onlinelibrary.wiley.com/doi/full/10.1096/fj.11-190868)
- [Functional cross-talk between allosteric effects of activating and inhibiting ligands underlies PKM2 regulation](https://elifesciences.org/articles/45068)

M32K25 structural alphabet was sourced from: https://github.com/AllosterIt/M32K25/tree/master

### Objective

Apply structural alphabet and information-theoretical approach to quantify convergence across multiple conformational states in molecular dynamics trajectories for a given sequence.

### High-level methodology

The pipeline reduces a trajectory to a string of structural-alphabet letters, quantifies the statistical coupling between positions with normalized mutual information, and compares coupling matrices across trajectory blocks with a matrix overlap.

```
  ┌───────────────────────────┐
  │   MD trajectory (xyz)      │   frames × atoms × 3
  │   T frames of Cα atoms     │
  └─────────────┬─────────────┘
                │  sliding window of 4 consecutive Cα atoms
                │  Kabsch-aligned to each M32K25 reference (min RMSD)
                ▼
  ┌───────────────────────────┐
  │   Structural encoding      │   encode_traj()
  │   frames × SA-letter string│   e.g. "A G K ... Y"
  └─────────────┬─────────────┘
                │  split into N trajectory blocks
                │  per block: letter counts → p(a), p(b), p(a,b)
                ▼
  ┌───────────────────────────┐
  │   Normalized mutual info   │   compute_nMI_for_traj_blocks()
  │   N blocks of L × L nMI    │   Iᵢⱼ = Hᵢ + Hⱼ − Hᵢⱼ, bias-corrected
  └─────────────┬─────────────┘
                │  pairwise matrix square-root distance (Hess 2002)
                ▼
  ┌───────────────────────────┐
  │   Matrix overlap Ω(A,B)    │   compute_matrix_overlap()
  │   block-vs-block similarity│   Ω ∈ [0, 1], 1 = identical coupling
  └───────────────────────────┘
```

### Key theoretical elements

The three core computations are summarized below.

#### 1. Kabsch alignment

Each fragment of four consecutive Cα atoms is encoded by finding the M32K25 reference letter that minimizes the optimally-superposed RMSD. Given a reference fragment $P$ and a candidate fragment $Q$ (each $\in \mathbb{R}^{4 \times 3}$), both are first mean-centered,

$$
P_c = P - \bar{P}, \qquad Q_c = Q - \bar{Q},
$$

and the cross-covariance matrix is decomposed with the singular value decomposition,

$$
H = P_c^{\top} Q_c, \qquad H = U \Sigma V^{\top}.
$$

The optimal rotation that superposes $Q_c$ onto $P_c$ is computed,

$$
d = \left(\det\left(V U^{\top}\right)\right), \qquad
R = V \, diag(1, 1, d) \, U^{\top},
$$

and the fragment is assigned the letter minimizing the resulting root-mean-square deviation over $N = 4$ atoms,

$$
\mathrm{RMSD}(P, Q) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \left\lVert p_i^{\,c} - R\, q_i^{\,c} \right\rVert^2 }.
$$

#### 2. Normalized mutual information

For any two encoding columns $i$ and $j$ (positions along the sequence, sampled over trajectory frames), letter counts are converted to marginal distributions $p_i(a)$, $p_j(b)$ and the joint distribution $p_{ij}(a, b)$. The Shannon entropy (in bits) and mutual information are

$$
H(X) = -\sum_{a} p(a) \log_2 p(a), \qquad
I_{ij} = H_i + H_j - H_{ij}.
$$

A finite-size correction accounts for the bias from the number of observed states ($N_i$, $N_j$ marginal, $N_{ij}$ joint) given $M$ frames, and the result is normalized by the joint entropy,

$$
\varepsilon_{ij} = \frac{N_{ij} - N_i - N_j + 1}{2 M}, \qquad
I^{\mathrm{norm}}_{ij} = \frac{I_{ij} - \varepsilon_{ij}}{H_{ij}} \quad (H_{ij} \neq 0).
$$

#### 3. Matrix overlap

Pairwise similarity between two symmetric positive semi-definite nMI matrices $A$ and $B$ is measured following Hess, *Phys. Rev. E* (2002), as a normalized distance between their matrix square roots,

$$
\Omega(A, B) = 1 - \frac{\sqrt{tr\!\left[\left(A^{1/2} - B^{1/2}\right)^2\right]}}{\sqrt{tr(A) + tr(B)}},
$$

which yields a value in $[0, 1]$, where $1$ indicates identical matrices.

### Potential applications

Use of structural alphabets and information theoretic approaches allows the following:
- Capturing observed conformational states from molecular dynamics simulations using string-based (e.g. MK32K25 alphabet) encoding.
- Encoding can be used as dynamic-based fingerprints to characterise conformational states or as input to machine learning models, e.g, sequence-based transformers.
- Compare conformation states of un-bound and interacting proteins and identify residues that are crucial in dynamic or allosteric interactions.

### Installation

Initialize `conda` environment:

```bash
conda create -n saencode -c conda-forge python=3.11
conda activate saencode
```

From the `saencode` directory, run the following command to install the package with optional dependencies:

```bash
pip install -e ".[dev]"
```

Run tests from a local folder:

```bash
pytest
```

### Package structure

`SAEncode` package has the following structure:

```
saencode/
├── __init__.py            # Public API: re-exports the functions below
├── encode.py              # Trajectory → structural-alphabet encoding (Kabsch alignment)
├── letters.py             # M32K25 alphabet: loads the 25 reference fragments at import
├── mutual_information.py   # Normalized mutual information between encoding positions
├── matrix_overlap.py      # Similarity (overlap) between two nMI matrices
└── alphabet/              # The 25 M32K25 reference fragments (A.pdb … Y.pdb)

tests/
└── test_saencode.py       # Validates encoding, nMI, and matrix overlap

tutorial/
└── example.ipynb          # End-to-end worked example 
```

Package implements and showcases fundamental aspects of:
- Applied linear algebra to compute structural comparisons and information-theoretic calculations.
- Packaging and deployment of a scientific software.
- Provided documentation through `README.md` and `numpy` docstring notation.
- Code reusability through examples and clear import structure.
- Ensuring code validity through testing.

### Tutorial

A high-level usage provided below:

```python
# Import libraries
import mdtraj as mdt
from MDAnalysisData import datasets
import saencode

# Load trajectory
traj = mdt.load("traj.xtc")

# Encode a trajectory
encoding = saencode.encode_traj(traj.xyz)

# Compute normalized mutual information for a given number of trajectory block
traj_nmi_blocks = saencode.compute_nMI_for_traj_blocks(encoding, number_of_blocks=5)

# Compute similarity (matrix overlap) between two normalized mutual information matrices
similarity = saencode.compute_matrix_overlap(traj_nmi_blocks)
```

A full visualized example can be found here: `tutorial/example.ipynb`

### Note on AI-assisted development

**2026 11 23:** Claude Opus 4.8 was employed to optimize the following:
- Review dependencies, setup `.toml`, and review imports.
- Improve `importlib` usage when loading M32K25 alphabet.
- Review docstrings and typing information.
- Update `README.md` file with initial Latex formulae.
- Update `README.md` file with project structure.
- Update `README.md` file with visual workflow.

### Future improvements:

Notes on improving package:
- Optimize computation of key computations and reduce dependencies on the looping.
- Implement tools to identify allosteric residue positions.
- Improve code testing.
