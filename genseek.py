import ase
""" The workflow"""

# Read the example geometry file or read SMILES code


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
