#!/usr/bin/env python

import numpy as np
import os.path
import sys

from scipy.spatial.distance import cdist


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

# Radii of back-back and base-base bounds
R_BACK       = 0.8
R_BASE       = 0.8

conf_counter = 0

[n_nucl, n_strands] = [int(x) for x in file_top.readline().split()]

idx_nucl  = 0
s_to_nucl = [[] for i in range(n_strands)]

nucl_line = file_top.readline()
timeline  = file_traj.readline()

while nucl_line:
	idx_strand = int(nucl_line.split()[0]) - 1
	s_to_nucl[idx_strand].append(idx_nucl)
	
	idx_nucl  += 1
	nucl_line  = file_top.readline()

while timeline:
	conf_counter += 1
	if conf_counter % 100 == 0: print("\033[1;34mProcessed %d configurations\033[0m" % conf_counter)

	box          = [float(x) for x in file_traj.readline().split()[2:]]
	[Et, Ep, Ek] = [float(x) for x in file_traj.readline().split()[2:5]]
	
	backs        = []
	bases        = []

	for idx_nucl in range(n_nucl):
		nucl_line = file_traj.readline().split()

		ci   = [float(x) for x in nucl_line[0:3]]
		a1   = [float(x) for x in nucl_line[3:6]]
		a3   = [float(x) for x in nucl_line[6:9]]

		ci   = np.asarray(ci)
		a1   = np.asarray(a1)
		a3   = np.asarray(a3)
		a2   = np.cross(a3,a1)

		base = ci + a1*POS_BASE
		back = ci + a1*POS_MM_BACK1 + a2*POS_MM_BACK2
		#back = ci + a1*POS_BACK
		
		bases.append(base)
		backs.append(back)

	bases = np.asarray(bases)
	backs = np.asarray(backs)

	s_list = []

	# Locate pairs of bound strands
	for idx_s1, s1 in enumerate(s_to_nucl[:-1]):
		for idx_s2, s2 in enumerate(s_to_nucl[idx_s1+1:]):
			backs1 = backs[s1,:]
			bases1 = bases[s1,:]

			backs2 = backs[s2,:]
			bases2 = bases[s2,:]

			b_back = np.min(cdist(backs1, backs2)) < R_BACK
			b_base = np.min(cdist(bases1, bases2)) < R_BASE
			
			if (b_back | b_base): s_list.append([idx_s1, idx_s1+idx_s2+1])

	s_list.extend([[idx_s] for idx_s in range(n_strands)])

	# Build fragment list from pairs involving common strands
	frags  = []

	while len(s_list) > 0:
		first, rest = s_list[0], s_list[1:]
		first       = set(first)
		lf          = -1
		
		while len(first) > lf:
			lf    = len(first)
			rest2 = []
			
			for r in rest:
				if len(first.intersection(set(r))) > 0: first |= set(r)
				else: rest2.append(r)
			
			rest = rest2
		
		frags.append(list(first))
		s_list = rest

	# Discard nucleotides from smaller fragments
	nucls       = [[idx_nucl for s in f for idx_nucl in s_to_nucl[s]] for f in frags]
	nucls_main  = max(nucls, key=len)

	backs       = backs[nucls_main, :]

	# Translate center-of-mass back to the origin
	center      = np.mean(backs, axis=0)

	backs      -= center
	points      = backs.T
	
	# Perform Principal Component Analysis of the configuration by singular-value decomposition
	P, D, Q     = np.linalg.svd(points)
	
	# P is the rotation matrix expressing the covariance matrix eigenvectors in the reference frame
	rot         = np.roll(P, 2, axis=1)
	inv_rot     = rot.T
	
	data_output = points
	
	# Rotate/save configuration so that the covariance-aligned OBB matches the frame with y,z the minimum,maximum dispersion axes
	data_output = np.dot(inv_rot, data_output)
	
	for i, nucl in enumerate(backs):
		file_output.write(" ".join(str(x) for x in data_output[:,i]))
		file_output.write("\n")

	file_output.write("\n\n")

	timeline = file_traj.readline()

file_traj.close()
file_top.close()
file_output.close()

print("\033[1;32mOutput printed to '%s'\033[0m" % path_output)
