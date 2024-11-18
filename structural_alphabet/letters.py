# Import libraries
import mdtraj as mdt

# Define a set of letters
LETTERS = [
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y'
]

# Load letter CA XYZ in a dictionary
M32K25 = {}
for letter in LETTERS:
    letter_pdb = mdt.load_pdb(filename="./alphabet/{}.pdb".format(letter))
    M32K25[letter] = letter_pdb.xyz[0]
