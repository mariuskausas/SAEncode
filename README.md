# structural_alphabet
Structural encoding and mutual information analysis of molecular dynamics trajectories.

Python implementation of "Detection of allosteric signal transmission by information-theoretic analysis of protein dynamics" by Alessandro Pandini et al, 2012:
- Paper link: https://faseb.onlinelibrary.wiley.com/doi/full/10.1096/fj.11-190868
- Structural Alphabet was sourced from: https://github.com/AllosterIt/M32K25/tree/master

## Required  third-party packages

```
mdtraj
numpy
scipy
tqdm
```

## Installation

```bash
# In the top level of the package
pip install .
```

## Tutorial

```python
# Import libraries
import structural_alphabet
```

```python
# Fetch an example molecular dynamics trajectory
from MDAnalysisData import datasets
adk = datasets.fetch_adk_equilibrium()
```

```python
# Pre-process molecular dynamics trajectory and extract only CA
top = mdt.load_psf(adk["topology"])
ca_atoms = top.select("name CA")
traj = mdt.load_dcd(
    filename=adk["trajectory"], 
    top=adk["topology"], 
    atom_indices=ca_atoms
    )
```

```
# More about the mdtraj trajectory object
<mdtraj.Trajectory with 4187 frames, 214 atoms, 214 residues, and unitcells at 0x7f31792ecc80>
```

```python
# Encode a trajectory
encoding = structural_alphabet.encode_traj(traj.xyz)
```

```
# traj_encoding
array([['B', 'G', 'A', ..., 'U', 'U', 'W'],
       ['E', 'G', 'G', ..., 'U', 'U', 'U'],
       ['E', 'D', 'I', ..., 'U', 'U', 'U'],
       ...,
       ['E', 'G', 'G', ..., 'U', 'U', 'U'],
       ['E', 'G', 'G', ..., 'U', 'U', 'W'],
       ['B', 'G', 'L', ..., 'U', 'U', 'W']], dtype='<U1')
```

```python
# Compute normalized mutual information for a given trajectory block
traj_nmi = structural_alphabet.compute_nMI_for_traj_block(encoding=encoding)
```

```python
# Compute similarity (matrix overlap) between two normalized mutual information matrices
similarity = structural_alphabet.compute_matrix_overlap([traj_nmi_1, traj_nmi_2])
```
