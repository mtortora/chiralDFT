#!/usr/bin/env python

import numpy as np
import os.path
import sys


# Input/output
if   len(sys.argv) == 3: path_output = os.path.splitext(sys.argv[1])[0] + ".in"
elif len(sys.argv) == 4: path_output = sys.argv[3]
else:
	print("\033[1;31mUsage is %s trajectory topology [output]\033[0m" % sys.argv[0])
	sys.exit()

path_traj    = sys.argv[1]
path_top     = sys.argv[2]

file_traj    = open(path_traj, mode="r")
file_top     = open(path_top, mode="r")
file_output  = open(path_output, mode="w")

POS_BACK     = -0.4
POS_BASE     = 0.4
POS_MM_BACK1 = -0.3400
POS_MM_BACK2 = 0.3408

conf_counter = 0

file_top.readline()

nucleotides  = file_top.readlines()
timeline     = file_traj.readline()

while timeline:
	conf_counter += 1
	if conf_counter % 100 == 0: print("\033[1;34mProcessed %d configurations\033[0m" % conf_counter)

	box          = [float(x) for x in file_traj.readline().split()[2:]]
	[Et, Ep, Ek] = [float(x) for x in file_traj.readline().split()[2:5]]
	
	cm_pos       = []

	for nucleotide in nucleotides:
		nucleo_line = file_traj.readline().split()

		ci   = [float(x) for x in nucleo_line[0:3]]
		a1   = [float(x) for x in nucleo_line[3:6]]
		a3   = [float(x) for x in nucleo_line[6:9]]

		ci   = np.asarray(ci)
		a1   = np.asarray(a1)
		a3   = np.asarray(a3)
		a2   = np.cross(a3,a1)

		base = ci + a1*POS_BASE
		back = ci + a1*POS_MM_BACK1 + a2*POS_MM_BACK2
		#back = ci + a1*POS_BACK
		
		cm_pos.append(back)

	# Translate center-of-mass back to the origin
	cm_pos      = np.asarray(cm_pos)
	center      = np.mean(cm_pos, axis=0)

	cm_pos     -= center
	points      = cm_pos.T

	# Perform Principal Component Analysis of the configuration by singular-value decomposition
	P, D, Q     = np.linalg.svd(points)

	# P is the rotation matrix expressing the covariance matrix eigenvectors in the reference frame
	rot         = np.roll(P, 2, axis=1)
	inv_rot     = rot.T

	data_output = points

	# Rotate/save configuration so that the covariance-aligned OBB matches the frame with y,z the minimum,maximum dispersion axes
	data_output = np.dot(inv_rot, data_output)
		
	for i,nucleotide in enumerate(nucleotides):
		file_output.write(" ".join(str(x) for x in data_output[:,i]))
		file_output.write("\n")

	file_output.write("\n\n")

	timeline = file_traj.readline()

file_traj.close()
file_top.close()
file_output.close()

print("\033[1;32mOutput printed to '%s'\033[0m" % path_output)