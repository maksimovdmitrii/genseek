from ase.io import read, write
from optparse import OptionParser

from ase import neighborlist
from ase.build import molecule
from scipy import sparse

from genseek.modules import detect_rotatble

parser = OptionParser()
parser.add_option("-g", "--geometry", dest="geometry",
                  help="Geometry")
parser.add_option("-f", "--format", dest="format",
                  help="Format of the geometry")
(options, args) = parser.parse_args()


""" The workflow"""

# Read the example geometry file or read SMILES code
atoms = read(options.geometry, format=options.format)
print(atoms)

# Preface
## Identify the connectivity of the geometry
cutOff = neighborlist.natural_cutoffs(atoms)
neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=False)
neighborList.update(atoms)
matrix = neighborList.get_connectivity_matrix()
list_of_torsions = detect_rotatble(matrix)

print(list_of_torsions)
## Identify the obstacles (Fixed frame)
## Specify number of molecules 
## Identify the periodic boundary conditions (PBC)

# Start generating of the ensembles
## Check for intermolecular clashes
## Check for intramolecular clashes
## Check for PBC clashes

# Optional Preexploration of the conformational space
## RMSD blacklisting
## Internal degrees of freedom blacklisting
## SOAP blacklisting

# Potential evaluation
## Blacklist check
## Run Minimization
## Blacklist check

# Next Trial
