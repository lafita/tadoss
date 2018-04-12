import argparse as ap
import warnings

import Bio.PDB as bp
from Bio import BiopythonWarning


# All the argument parsing logic
parser = ap.ArgumentParser()

parser.add_argument("input", help="path to the input structure file", type=str)
parser.add_argument("output", help="name of the output structure file", type=str)

parser.add_argument("-w", "--warn", help="print warning messages",
                    action='store_true')

args = parser.parse_args()

# Silent BioPython warnings
if not args.warn:
    warnings.simplefilter('ignore', BiopythonWarning)

# Parse the PDB file into a structure object
parser = bp.PDBParser()
structure = parser.get_structure("input", args.input)

# Only select amino acid residues (protein) to be written in the output
class ProteinSelect(bp.Select):
    def accept_residue(self, residue):
        if residue.has_id('CA') :
            return 1
        else:
            return 0

# Write the output file
io = bp.PDBIO()
io.set_structure(structure)
io.save(args.output, ProteinSelect())
