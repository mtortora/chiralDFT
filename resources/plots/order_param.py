#!/usr/bin/env python3

import sys
import os.path

import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv) != 5:
	print("\033[1;31mUsage is %s order_file l mp m\033[0m" % sys.argv[0])
	sys.exit()


# Input/output
path_in = sys.argv[1]

l       = int(sys.argv[2])
mp      = int(sys.argv[3])
m       = int(sys.argv[4])

if not os.path.isfile(path_in):
	print("\033[1;31mCouldn't open file %s - aborting\033[0m" % path_in)
	sys.exit()

if (mp < -l) | (mp > l) | (m < -l) | (m > l):
	print("\033[1;31mInconsistent set of indices (must have -l <= mp,m <= l)\033[0m")
	sys.exit()


# Format c++ complex output
data     = np.genfromtxt(path_in, dtype=str)
data     = np.vectorize(lambda x: complex(*eval(x)))(data)

idx_l    = 2*l * (l+1)*(2*l+1)/3 - l*(2*l+1)
idx_mp   = (2*l+1)*(mp+l)
idx_m    = m+l

idx      = idx_l + idx_mp + idx_m + 1

# Plot
eta_grid = np.real(data[:,0])

real     = np.real(data[:,int(idx)])
imag     = np.imag(data[:,int(idx)])

plt.plot(eta_grid, real, label=r'$\Re(D^%d_{%d,%d})$' % (l,mp,m))
plt.plot(eta_grid, imag, label=r'$\Im(D^%d_{%d,%d})$' % (l,mp,m))

plt.legend(loc=4)

plt.show()
