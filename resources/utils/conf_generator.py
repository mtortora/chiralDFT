#!/usr/bin/env python

import numpy as np
import os.path
import sys


if len(sys.argv) != 3:
	print("\033[1;31mUsage is %s wireframe num_box\033[0m" % sys.argv[0])
	sys.exit()

path_wire    = sys.argv[1]
num_box      = sys.argv[2]

path_ref     = os.path.splitext(path_wire)[0] + ".ref"
path_cnf     = os.path.splitext(path_wire)[0] + ".cnf"
path_ind     = os.path.splitext(path_wire)[0] + ".ind"

data_wire    = np.genfromtxt(path_wire)
num_box      = int(num_box)

# Prune duplicate vertices
data_wire    = np.vstack({tuple(vertex) for vertex in data_wire})
data_wire    = data_wire.T

# Work out semi-oriented bounding box parameters
norms        = np.linalg.norm(data_wire, axis=0)
r_max        = np.max(norms)
r_integ      = 2 * r_max

radial_norms = np.linalg.norm(data_wire[:2,:], axis=0)
l_rc         = np.max(radial_norms)

l_z          = np.max(data_wire[2,:]) - np.min(data_wire[2,:])
l_zh         = l_z / 2.
l_bnd_zh     = l_zh / num_box

scale        = l_zh - l_bnd_zh
rescaled     = np.linspace(-scale, scale, num_box)

xaxis_ref    = np.array([1.,0.,0.])
zaxis_ref    = np.array([0.,0.,1.])

centers_ref  = np.outer(rescaled, zaxis_ref)

step_counter = 0
box_overlaps = np.zeros([num_box, num_box])


# Rotation matrix from axis and angle (Euler-Rodrigues)
def rotation_matrix(axis, angle):
	axis    = np.asarray(axis)
	axis   /= np.linalg.norm(axis)
	
	a       = np.cos(angle/2.)
	b, c, d = -axis*np.sin(angle/2.)
	
	aa, bb, cc, dd         = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	
	return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
         	         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
         	         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


# Semi-oriented bounding box overlap test
def overlap_box(r_sep_box, axis1, axis2):
	dim_box = np.array([l_rc, l_rc, l_bnd_zh])
	
	cos_t   = np.dot(axis1, axis2)
	sin_t   = np.sqrt(1 - cos_t**2)
	
	box1_x  = np.cross(axis1, axis2) / sin_t
	box1_y  = np.cross(axis1, box1_x)
	box1    = np.array([box1_x, box1_y, axis1])
	
	theta   = np.arccos(cos_t)
	r_proj  = np.dot(box1, r_sep_box)
	
	rot     = rotation_matrix(xaxis_ref, theta)
	abs_rot = np.abs(rot)
	
	for i in range(3):
		r1 = dim_box[i]
		r2 = np.dot(dim_box, abs_rot[i,:])
		
		if np.abs(r_proj[i]) > r1 + r2: return False
	
	for i in range(1,3):
		r1 = np.dot(dim_box, abs_rot[:,i])
		r2 = dim_box[i]
		
		if np.abs(np.dot(r_proj, rot[:,i])) > r1 + r2: return False
	
	return True


# Generate configurations
while not np.any(box_overlaps):
	step_counter += 1
	
	angle         = np.random.uniform(0, 2*np.pi)
	r_sep         = (np.random.random(3)-0.5) * 2*r_integ
	
	rot_cnf       = rotation_matrix(xaxis_ref, angle)
	zaxis_cnf     = rot_cnf[:,2]
	
	centers_cnf   = np.outer(rescaled, zaxis_cnf)
	
	for idx_box_ref in range(num_box):
		for idx_box_cnf in range(num_box):
			vec_sep = r_sep + centers_cnf[idx_box_cnf,:] - centers_ref[idx_box_ref,:]
			box_overlaps[idx_box_ref, idx_box_cnf] = overlap_box(vec_sep, zaxis_ref, zaxis_cnf)

# Save randomised configuration in non-rotated form for easier box computations
wire_ref  = data_wire.copy()
wire_cnf  = data_wire.copy() + r_sep[:,np.newaxis]

# Recover orientation as the head vector of the configuration file
zaxis_cnf = zaxis_cnf[:,np.newaxis]
data_cnf  = np.hstack((zaxis_cnf, wire_cnf))

np.savetxt(path_ref, wire_ref.T)
np.savetxt(path_cnf, data_cnf.T)

np.savetxt(path_ind, box_overlaps)

print("\033[1;34mBoxes overlapped after %d steps\033[0m" % step_counter)

print("\033[1;32mReference configuration printed to '%s'\033[0m" % path_ref)
print("\033[1;32mRandomised configuration printed to '%s'\033[0m" % path_cnf)
print("\033[1;32mOverlap indices printed to '%s'\033[0m" % path_ind)
