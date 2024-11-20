# Import libraries
import os
import importlib.resources 
from pathlib import Path 
import mdtraj as mdt


# Define a set of letters
LETTERS = [
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y'
]

# Load letter CA XYZ in a dictionary
M32K25 = {}
with importlib.resources.path('saencode.alphabet', '') as resource_path: 
    data_path = Path(resource_path)
    path_to_pdbs = list(data_path.iterdir())
    for path_to_pdb in path_to_pdbs:
        file_basename = os.path.basename(path_to_pdb)
        if "pdb" in file_basename:
            letter = file_basename[0]
            letter_pdb = mdt.load_pdb(filename=path_to_pdb)
            M32K25[letter] = letter_pdb.xyz[0]
