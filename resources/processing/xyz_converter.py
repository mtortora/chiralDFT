#!/usr/bin/env python

import sys
import os.path


# Input/output
if len(sys.argv) in [2,3]:
	path_in  = os.path.realpath(sys.argv[1])

	if len(sys.argv) == 3: path_out = sys.argv[2]
	else: path_out = os.path.dirname(path_in) + "/trajectory.in"
	
else:
	print("\033[1;31mUsage is %s trajectory [output]\033[0m" % sys.argv[0])
	sys.exit()

file_in  = open(path_in,  mode="r")
file_out = open(path_out, mode="w")


# Parse trajectory
c_count  = 0
s_line   = file_in.readline()

while s_line:
	file_in.readline()

	c_count += 1
	n_s      = int(s_line.split()[0])

	for idx_s in range(n_s):
		xs = [float(x) for x in file_in.readline().split()[1:]]
		
		file_out.write(" ".join(str(x) for x in xs))
		file_out.write("\n")

	file_out.write("\n\n")

	s_line = file_in.readline()


file_out.close()

print("\033[1;32mPrinted %d configurations to '%s'\033[0m" % (c_count, path_out))
