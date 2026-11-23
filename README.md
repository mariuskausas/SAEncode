## Structural encoding and mutual information analysis of molecular dynamics trajectories.

### Reference

Python implementation of structural encoding using M32K25 and mutual information analysis of molecular dynamics trajectories based on the following references:

- [Detection of allosteric signal transmission by information-theoretic analysis of protein dynamics](https://faseb.onlinelibrary.wiley.com/doi/full/10.1096/fj.11-190868)
- [Functional cross-talk between allosteric effects of activating and inhibiting ligands underlies PKM2 regulation](https://elifesciences.org/articles/45068)

M32K25 structural alphabet was sourced from: https://github.com/AllosterIt/M32K25/tree/master

### Objective



### Approach

The pipeline reduces a trajectory to a string of structural-alphabet letters, quantifies the statistical coupling between positions with normalized mutual information, and compares coupling matrices across trajectory blocks with a matrix overlap. The three core computations are summarized below.

#### 1. Kabsch alignment

Each fragment of four consecutive Cα atoms is encoded by finding the M32K25 reference letter that minimizes the optimally-superposed RMSD. Given a reference fragment $P$ and a candidate fragment $Q$ (each $\in \mathbb{R}^{4 \times 3}$), both are first mean-centered,

$$
P_c = P - \bar{P}, \qquad Q_c = Q - \bar{Q},
$$

and the cross-covariance matrix is decomposed with the singular value decomposition,

$$
H = P_c^{\top} Q_c, \qquad H = U \Sigma V^{\top}.
$$

The optimal rotation that superposes $Q_c$ onto $P_c$ is recovered while correcting for an improper reflection,

$$
d = \operatorname{sign}\!\left(\det\!\left(V U^{\top}\right)\right), \qquad
R = V \, \operatorname{diag}(1, 1, d) \, U^{\top},
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
\Omega(A, B) = 1 - \frac{\sqrt{\operatorname{tr}\!\left[\left(A^{1/2} - B^{1/2}\right)^2\right]}}{\sqrt{\operatorname{tr}(A) + \operatorname{tr}(B)}},
$$

which yields a value in $[0, 1]$, where $1$ indicates identical matrices.


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

### Note on AI-drive development

**2026 11 23:** Claude Opus 4.8 was employed to optimize the following
- Review dependencies, setup `.toml`, and review imports.
- Improve `importlib` usage when loading M32K25 alphabet.
- Review docstrings and typing information.
- Update `README.md` file with initial Latex formulaes.
- Update `README.md` file with project structure.
