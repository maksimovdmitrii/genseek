from ase.io import read, write
from optparse import OptionParser

from genseek.modules import *
import numpy as np
import sys

from random import randint

parser = OptionParser()
parser.add_option("-g", "--geometry", dest="geometry",
                  help="Geometry")
parser.add_option("-f", "--geometryformat", dest="format",
                  help="Format of the geometry")
parser.add_option("--fixedgeometry", dest="fixedgeometry",
                  help="Format of the fixed geometry")
parser.add_option("--fixedgeometryformat", dest="fixedgeometryformat",
                  help="Format of the fixed geometry")
(options, args) = parser.parse_args()


""" The workflow"""

# Read the example geometry file or read SMILES code
atoms = read(options.geometry, format=options.format)
# print(atoms)

# Preface
## Identify the connectivity of the geometry
connectivity_matrix_full = create_connectivity_matrix(atoms, bothways=True)
connectivity_matrix_isolated = create_connectivity_matrix(atoms, bothways=False)
list_of_torsions = detect_rotatble(connectivity_matrix_isolated)
# print(list_of_torsions)

## Identify the obstacles (Fixed frame) and PBC
# Minimum Image Convention
mic = False
if options.fixedgeometry is not None:
	fixed_frame = read(options.fixedgeometry, format=options.fixedgeometryformat)
	pbc = fixed_frame.get_cell()[:]
	if not np.array_equal(pbc, np.zeros([3, 3])):
		mic = True
		atoms.set_cell(pbc)
		atoms.set_pbc(True)

# Set atoms template to initial configuration
set_centre_of_mass(atoms, np.array([0, 0, 0]))
align_to_axes(atoms,  0, len(atoms)-1)
# Check for clashes

# connectivity_matrix = create_connectivity_matrix(atoms)
# print(connectivity_matrix)
# if not internal_clashes(atoms, connectivity_matrix_full):
# 	sys.exit('The Cell is too small! Increase the cell size!')

## Specify number of molecules 
molecules = [atoms.copy() for i in range(3)]

# Blacklist
blacklist = create_blacklist(molecules, list_of_torsions)
	
Trials = 10000
Trial = 0

while Trial<Trials:
	Trial+=1
	print("Start the new Trial ", Trial)
	vector = create_internal_vector(molecules, list_of_torsions)
	if not_in_blacklist(vector, blacklist):
		blacklist = add_to_blacklist(vector, blacklist)
		apply_to_molecules(molecules, vector, list_of_torsions, connectivity_matrix_isolated)
		# print(measure_molecules(molecules, list_of_torsions))
		if all_right(molecules, fixed_frame, connectivity_matrix_full, mic=mic):
			ensemble = atoms.copy()
			del ensemble[[atom.index for atom in atoms]]
			for molecule in molecules:
				ensemble+=molecule
			ensemble+=fixed_frame
			write("initial.in", ensemble, format="aims")
			break
		else:
			ensemble = atoms.copy()
			del ensemble[[atom.index for atom in atoms]]
			for molecule in molecules:
				ensemble+=molecule
			ensemble+=fixed_frame
			write("WithClash.in", ensemble, format="aims")
			print("Next trial, produced with clash")
			continue
	else:
		print("Next trial, found in blacklist")
		continue

# Write the enesemble into file 
# add the fixed frame also
# Extend method can be used instead this:






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
