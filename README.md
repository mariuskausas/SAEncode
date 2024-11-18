# structural_alphabet
Python implementation of Structural Alphabet encoding of molecular dynamics trajectories.

## To-do list:
- Test the tutorial.
- Provide required installation documentation.
- Provide a theoretical background.
- Showcase how to perform a simple encoding and mutual inforamtion analysis.
- Write docstrings.

## Structural Alphabet

Alphabet letter coordinates downloaded from: https://github.com/AllosterIt/M32K25/tree/master

## Tutorial

```python
# Import libraries
import mdtraj as mdt
from MDAnalysisData import datasets

from .letters import M32K25
from .encode import encode_traj
from .mutual_information import compute_nMI_for_traj_block
```

```python
# Download data
adk = datasets.fetch_adk_equilibrium()
```

```python
# Load dataset
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
traj_encoding = encode_traj(traj_xyz=traj.xyz)
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
traj_nMI = compute_nMI_for_traj_block(encoding=traj_encoding)
```
