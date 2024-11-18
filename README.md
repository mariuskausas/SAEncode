# structural_alphabet
Python implementation of Structural Alphabet encoding of molecular dynamics trajectories.

## To-do list:
- Provide required installation documentation.
- Provide a theoretical background.
- Showcase how to perform a simple encoding and mutual inforamtion analysis.
- Write docstrings.

## Structural Alphabet

Structural Alphabet is source from: https://github.com/AllosterIt/M32K25/tree/master

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
# Compute normalized Mutual Information
traj_nmi = structural_alphabet.compute_nMI_for_traj_block(encoding=encoding)
```
