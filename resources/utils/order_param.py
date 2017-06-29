#!/usr/bin/env python3

import sys
import os.path

import numpy as np


if len(sys.argv) != 3:
	print("\033[1;31mUsage is %s order_file l_max\033[0m" % sys.argv[0])
	sys.exit()


# Input/output
path_in  = sys.argv[1]
path_out = os.path.splitext(path_in)[0]

l_max    = int(sys.argv[2])
size     = (l_max+1)*(2*l_max+1) + 2*l_max * (l_max+1)*(2*l_max+1)/3

if not os.path.isfile(path_in):
	print("\033[1;31mCouldn't open file %s - aborting\033[0m" % path_in)
	sys.exit()


# Format c++ complex output
data_in  = np.genfromtxt(path_in, dtype=str)
data_in  = np.vectorize(lambda x: complex(*eval(x)))(data_in)

eta_grid = np.real(data_in[:,0])

for l in range(l_max+1):
	for mp in range(-l,l+1):
		for m in range(-l,l+1):

			idx_l    = 2*l * (l+1)*(2*l+1)/3 - l*(2*l+1)
			idx_mp   = (2*l+1)*(mp+l)
			idx_m    = m+l

			idx      = idx_l + idx_mp + idx_m + 1

			real     = np.real(data_in[:,int(idx)])
			imag     = np.imag(data_in[:,int(idx)])

			data_out = np.array([eta_grid, real, imag]).T
			file_out = "%s_D%d%d%d.res" % (path_out,l,mp,m)

			np.savetxt(file_out, data_out)

print("\033[1;32mPrinted %d order parameters to '%s_Dxxx'\033[0m" % (size, path_out))
