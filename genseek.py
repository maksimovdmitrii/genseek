""" Generate and Search"""

__author__ = "GenSeek"
__copyright__ = "Copyright (C) 2020 Dmitrii Maksimov"
__license__ = "Public Domain"
__version__ = "0.1.1"

from optparse import OptionParser

from genseek.blacklist import *
from genseek.outputs import *
from genseek.structure import *

import numpy as np
import sys
from random import randint, random, uniform

parser = OptionParser()
parser.add_option("-t", "--test")
(options, args) = parser.parse_args()

""" Prefase"""
dirs = Directories()
output = Output("report.out")
parameters = load_parameters("parameters.json")
workflow = Workflow()
structure = Structure(parameters)
fixed_frame = Fixed_frame(parameters)
blacklist = Blacklist(structure)


while workflow.trials < parameters["trials"]:
	workflow.trials += 1
	print("New Trial", workflow.trials)
	output.write("Start the new Trial {}\n".format(workflow.trials))
	# Generate the vector in internal degrees of freedom
	configuration = structure.create_configuration(parameters)
	if blacklist.not_in_blacklist(configuration):
		blacklist.add_to_blacklist(configuration)
		structure.apply_configuration(configuration)
		if all_right(structure, fixed_frame):
			dirs.create_directory()
			ensemble = merge_together(structure, fixed_frame)
			dirs.save_to_directory(ensemble, parameters)
	else:
		output.write("Next trial, found in blacklist")
		continue

# Write the enesemble into file 
# add the fixed frame also
# Extend method can be used instead this:





## Identify the periodic boundary conditions (PBC)

# Start generating of the ensembles
## Check for intermolecular clashes
## Check for intramolecular clashes
## Check for PBC clashes
# 
# Optional Preexploration of the conformational space
## RMSD blacklisting
## Internal degrees of freedom blacklisting
## SOAP blacklisting

# Potential evaluation
## Blacklist check
## Run Minimization
## Blacklist check

# Next Trial
