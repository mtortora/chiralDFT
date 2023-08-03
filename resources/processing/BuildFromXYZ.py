#!/usr/bin/env python

import sys
import os.path


# Input/output
if len(sys.argv) in [2,3]:
	path_in = os.path.realpath(sys.argv[1])

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
		linedata = file_in.readline().split()
		
		atype = linedata[0]
		charge = float(linedata[-1])

		xs = [x for x in linedata[1:]]
		
		if atype == "CH1":
			file_out.write("0 ")
		elif atype == "CH2":
			file_out.write("1 ")
		elif atype == "CH3":
			file_out.write("2 ")
		elif atype == "CH2r":
			file_out.write("3 ")
		elif atype == "CR1":
			file_out.write("4 ")
		elif atype == "N":
			file_out.write("5 ")
		elif atype == "NR":
			file_out.write("6 ")
		elif atype == "NT":
			file_out.write("7 ")
		elif atype == "O":
			file_out.write("8 ")
		elif atype == "OA":
			file_out.write("9 ")
		elif atype == "OM":
			file_out.write("10 ")
		elif atype == "S":
			file_out.write("11 ")
		elif (atype == "H") & (charge != 0.):
			file_out.write("12 ")
		
		if (atype != "H") | (charge != 0.):
			file_out.write(" ".join(x for x in xs))
			file_out.write("\n")

	file_out.write("\n\n")

	s_line = file_in.readline()

file_in.close()
file_out.close()

print("\033[1;32mPrinted %d configurations to '%s'\033[0m" % (c_count, path_out))
