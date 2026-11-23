# Import libraries
import importlib
from importlib.resources import files
import mdtraj as mdt


# Define a set of letters
LETTERS = [
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y'
]

# Load letter CA XYZ in a dictionary
M32K25 = {}
data_path = files('saencode.alphabet')
for p in data_path.iterdir():
    if p.name.endswith('.pdb'):
        letter = p.name[0]
        with importlib.resources.as_file(p) as fp:
            M32K25[letter] = mdt.load_pdb(filename=str(fp)).xyz[0]
