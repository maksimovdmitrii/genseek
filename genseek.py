from ase.io import read, write
from optparse import OptionParser




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
