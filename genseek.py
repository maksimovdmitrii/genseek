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
print(atoms)

# Preface
## Identify the connectivity of the geometry
list_of_torsions = detect_rotatble(atoms)
print(list_of_torsions)

# Set atoms template to initial configuration
set_centre_of_mass(atoms, np.array([0, 0, 0]))
# for torsion in list_of_torsions:
# 	fixed_indices = carried_atoms(atoms, torsion)
# 	atoms.set_dihedral(angle=180,
# 							  a1=torsion[0],
# 							  a2=torsion[1],
# 							  a3=torsion[2],
# 							  a4=torsion[3],
# 							  indices=fixed_indices)
# print(atoms.get_moments_of_inertia(vectors=True))
print(measure_quaternion(atoms, 0, len(atoms)-1))
align_to_axes(atoms,  0, len(atoms)-1)
print(measure_quaternion(atoms, 0, len(atoms)-1))
quaternion_set(atoms, produce_quaternion(0, np.array([1, 0, 0])), 0, len(atoms)-1)
print(measure_quaternion(atoms, 0, len(atoms)-1))

## Specify number of molecules 
molecules = [atoms.copy() for i in range(5)]
print(molecules)

for i in range(len(molecules)):
	# Set center of mass
	value = np.array([randint(-10,10),
					  randint(-10,10),
					  randint(-10,10)])
	set_centre_of_mass(molecules[i], value)
	# Set torsions
	for torsion in list_of_torsions:
		value = randint(0, 360)
		fixed_indices = carried_atoms(molecules[i], torsion)
		molecules[i].set_dihedral(angle=value,
								  a1=torsion[0],
								  a2=torsion[1],
								  a3=torsion[2],
								  a4=torsion[3],
								  indices=fixed_indices)
	# Set rotation
	quaternion_set(molecules[i], produce_quaternion(0, np.array([0, 0, 1])), 0, len(atoms)-1)


# Write the enesemble into file 
# add the fixed frame also
# Extend method can be used instead this:
ensemble = atoms.copy()
del ensemble[[atom.index for atom in ensemble]]
for molecule in molecules:
	ensemble+=molecule

## Identify the obstacles (Fixed frame)
if options.fixedgeometry is not None:
	fixed = read(options.fixedgeometry, format=options.fixedgeometryformat)
	ensemble+=fixed

write("initial.xyz", ensemble, format="xyz")


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
