#!/usr/bin/env python3

import sys
import os.path

import numpy as np


try:
	import spherical_functions as sf

except ImportError:
	print("\033[1;31mCould not find spherical_functions module - install via pip\033[0m")
	sys.exit()


if len(sys.argv) != 3:
	print("\033[1;31mUsage is %s biaxial_psi_file n_l\033[0m" % sys.argv[0])
	sys.exit()


# Input/output
path_in    = sys.argv[1]
n_l        = int(sys.argv[2])

name, ext  = os.path.splitext(path_in)
path_out   = name + "_w%d.res" % n_l

file_out   = open(path_out, mode="w")

if not os.path.isfile(path_in):
	print("\033[1;31mCouldn't open file %s - aborting\033[0m" % path_in)
	sys.exit()

data       = np.genfromtxt(path_in)

# Reshape
alphas     = data[:,0]
thetas     = data[:,1]
phis       = data[:,2]
psi        = data[:,3]

n_a0       = len(alphas[alphas==alphas[0]])
n_t0       = len(thetas[thetas==thetas[0]])
n_p0       = len(phis  [phis==phis[0]])

n_alpha    = int(np.sqrt(n_p0*n_t0/n_a0))
n_theta    = int(np.sqrt(n_p0*n_a0/n_t0))
n_phi      = int(np.sqrt(n_a0*n_t0/n_p0))

psi        = psi.reshape([n_alpha, n_theta, n_phi])

# Remove endpoints to avoid poles double-counting
alpha_grid = np.linspace(0,2.*np.pi, num=n_alpha, endpoint=False)
theta_grid = np.linspace(0,1.*np.pi, num=n_theta, endpoint=False)
phi_grid   = np.linspace(0,2.*np.pi, num=n_phi,   endpoint=False)

d_alpha    = 2.*np.pi / n_alpha
d_theta    = 1.*np.pi / n_theta
d_phi      = 2.*np.pi / n_phi

indices    = np.array([[l,mp,m] for l in range(n_l) for mp in range(-l,l+1) for m in range(-l,l+1)])
nw         = len(indices)

aves       = np.zeros(nw, dtype='complex128')

print("\033[1;36mIncluding %d Wigner coefficients\033[0m" % nw)


# Average order parameters
for i in range(n_alpha):
	for j in range(n_theta):
		for k in range(n_phi):
			alpha = alpha_grid[i]
			theta = theta_grid[j]
			phi   = phi_grid  [k]
			
			Ds    = sf.Wigner_D_element(phi, theta, alpha, indices)
			aves += np.conj(Ds) * psi[i,j,k] * np.sin(theta)*d_theta*d_phi*d_alpha

	print("\033[1;34mAveraged over %d out of %d iso-alpha surfaces\033[0m" % (i+1, n_alpha))

print("\033[1;32mCompleted averaging process\033[0m")

# Wigner decomposition
for i in range(n_alpha):
	for j in range(n_theta):
		for k in range(n_phi):
			alpha = alpha_grid[i]
			theta = theta_grid[j]
			phi   = phi_grid  [k]
			
			Ds    = sf.Wigner_D_element(phi, theta, alpha, indices)
			
			for idx in range(nw): Ds[idx] *= (2.*indices[idx][0]+1.)/(8*np.pi**2)
			
			psi_c = np.real((Ds*aves).sum())
			
			file_out.write(str(alpha) + " " + str(theta) + " " + str(phi) + " " + str(psi_c))
			file_out.write("\n")
		
		file_out.write("\n")

	print("\033[1;34mReconstructed %d out of %d iso-alpha surfaces\033[0m" % (i+1, n_alpha))

	file_out.write("\n")

file_out.close()

print("\033[1;32mOutput printed to '%s'\033[0m" % path_out)
