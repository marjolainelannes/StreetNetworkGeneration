#!/usr/bin/env python

import os
import sys
import platform

# Use this instead of os.system so that command
# failure terminates this script
# with corresponding return code.
def myrun(command):
	status = os.system(command)
	print(command)
	if status != 0:
		sys.exit(status)

####################################################
# Location of the work directory, on each
# machine and location of the final save
# directory on the RAID

# Import machine where the program is running
machine = platform.node()

#####################################################
# Read the list of processes
sys.argv.pop(0)
parameters=dict()
for arg in sys.argv:
	members=arg.split('=')
	parameters[members[0]]=members[1]
print(parameters)

# Extraction of the processes
N_cell = parameters["NCELL"]

# File destruction and creation
myfile = "road_width.py"
relative_path = N_cell

# lauching of the program
myrun("python " + myfile + " " + relative_path)
