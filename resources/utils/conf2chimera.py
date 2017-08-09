#!/usr/bin/env python

import numpy as np
import math as m
import os.path
import sys


# Load input files
if len(sys.argv) in [2,3,4]:
	display_boxes = False
	
	path_ref      = os.path.realpath(sys.argv[1])
	wire_ref      = np.genfromtxt(path_ref)
	
	path_com      = os.path.dirname(path_ref) + '/chimera.com'
	file_com      = open(path_com, mode="w")
	
	filenames     = [path_ref]
	wires         = [wire_ref]
	
	if len(sys.argv) > 2:
		path_cnf = sys.argv[2]
		data_cnf = np.genfromtxt(path_cnf)
		
		axis_cnf = data_cnf[0,:]
		wire_cnf = data_cnf[1:,:]
		
		if len(wire_cnf) != len(wire_ref):
			print("\033[1;31mInconsistent input files - aborting\033[0m")
			sys.exit()
		
		filenames.append(path_cnf)
		wires.append(wire_cnf)
		
		if len(sys.argv) == 4:
			display_boxes = True
			
			path_idx      = sys.argv[3]
			idx_ovrlp     = np.genfromtxt(path_idx)

else:
	print("\033[1;31mUsage is %s config_ref [config [indices]]\033[0m" % sys.argv[0])
	sys.exit()


# Global chimera presets
base_atom   = "H"

atom_color  = "0.2734,0.5078,0.7031"
atom_alpha  = 0.2

vdw_radius  = 0.5
ball_scale  = 1
atom_radius = vdw_radius*ball_scale

n_atoms     = len(wire_ref)
centers     = [np.mean(wires[idx_wire], axis=0) for idx_wire in range(len(wires))]

file_com.write("~bond\n")
file_com.write("represent b+s\n")
file_com.write("preset apply pub 3\n")

file_com.write("setattr m ballScale " + str(ball_scale) + "\n")
file_com.write("vdwdefine " + str(vdw_radius) + " @" + base_atom + "=\n")

file_com.write("colordef color_opaque " + atom_color + "\n")
file_com.write("colordef color_transp " + atom_color + "," + str(atom_alpha) + "\n")

if display_boxes: file_com.write("color color_transp " + base_atom + "\n")
else:             file_com.write("color color_opaque " + base_atom + "\n")


# Work out bounding parameters and relative orientations
center = (np.min(wire_ref, axis=0) + np.max(wire_ref, axis=0)) / 2.

z_min  = np.min(wire_ref[:,2])
z_max  = np.max(wire_ref[:,2])

if display_boxes:
	n_bnd     = len(idx_ovrlp)
	range_bnd = np.linspace(z_min, z_max, n_bnd+1)
	
	bnd_ovrlp = np.nonzero(idx_ovrlp)
	atoms_in  = [[] for idx_bnd in range(n_bnd)]
	
	# Assign atoms to z-stacked bounding boxes on reference configuration
	for idx,vertex in enumerate(wire_ref):
		for idx_bnd in range(n_bnd):
			vtx_in_bnd = (range_bnd[idx_bnd] <= vertex[2] <= range_bnd[idx_bnd+1])
			if vtx_in_bnd: atoms_in[idx_bnd].append(base_atom + str(idx+1))

	atoms_in = [",".join(atom for atom in atoms_in[idx_bnd]) for idx_bnd in range(n_bnd)]

# Compute spherocylinders only if box display is off
else:
	radial_norms = np.linalg.norm((wire_ref-center)[:,:2], axis=1)
	
	l_cr         = np.max(radial_norms)
	r_axis       = atom_radius + l_cr


# Process wireframes as collections of base_atoms
for idx_wire,wire in enumerate(wires):
	path_xyz = filenames[idx_wire] + '.xyz'
	file_xyz = open(path_xyz, mode="w")
	
	file_xyz.write(str(n_atoms) + "\n")
	file_xyz.write("mesogen " + str(idx_wire) + "\n")
	
	file_com.write("center #" + str(idx_wire) + "\n")
	
	# Write xyz file
	for vertex in wire:
		vertex = list(vertex)
		vertex.insert(0, base_atom)
		
		file_xyz.write(" ".join(str(x) for x in vertex))
		file_xyz.write("\n")

	if display_boxes:
		for idx_bnd in range(n_bnd):
			idx_box = idx_wire*n_bnd + idx_bnd + 2
			
			# Work out bounding boxes through local molecular maps
			file_com.write("molmap #")
			file_com.write(str(idx_wire) + ":@")
			file_com.write(atoms_in[idx_bnd])
			file_com.write(" " + str(atom_radius) + " modelId " + str(idx_box))
			file_com.write(" center ")
			file_com.write(",".join(str(x) for x in centers[idx_wire]))
			file_com.write(" edgePadding 0 showDialog 0 symmetry c4\n")
			
			file_com.write("volume #" + str(idx_box))
			file_com.write(" showOutlineBox 1 outlineBoxLinewidth 1.5 outlineBoxRgb")
			
			if idx_bnd in bnd_ovrlp[idx_wire]: file_com.write(" red\n")
			else:                              file_com.write(" black\n")

	file_xyz.close()
	
	print("\033[1;32mOutput %d printed to '%s'\033[0m" % tuple([idx_wire,path_xyz]))

if display_boxes:
	cur_dir = os.path.dirname(os.path.realpath(__file__))
	
	file_com.write("cd " + cur_dir + "\n")
	file_com.write("runscript prune_surfaces.py\n")


# Display cylindrical axes and spherical caps
else:
	for idx_wire,wire in enumerate(wires):
		file_com.write("shape cylinder radius ")
		file_com.write(str(r_axis))
		file_com.write(" center ")
		file_com.write(",".join(str(x) for x in center))
		file_com.write(" height " + str(z_max-z_min))
		file_com.write(" color color_transp")
		file_com.write("\n")
		
		for z in [z_min,z_max]:
			file_com.write("shape sphere radius ")
			file_com.write(str(r_axis))
			file_com.write(" center ")
			file_com.write(",".join(str(x) for x in center[:2]))
			file_com.write("," + str(z) + " color color_transp")
			file_com.write("\n")


# Rotate random configuration back
if len(wires) == 2:
	theta = np.arccos(axis_cnf[2])
	theta = m.copysign(theta, axis_cnf[1])
	
	file_com.write("turn x ")
	file_com.write(str(theta * 180/np.pi))
	file_com.write(" models #0")
	
	if display_boxes: file_com.write(",2-" + str(n_bnd+1))
	
	file_com.write(" center ")
	file_com.write(",".join(str(x) for x in centers[1]))
	file_com.write("\n")
	
	# Restore opacity of "overlapping" atoms
	if display_boxes:
		for idx_wire,wire in enumerate(wires):
			for idx_bnd in range(n_bnd):
				if idx_bnd in bnd_ovrlp[idx_wire]:
					file_com.write("transparency 0 #")
					file_com.write(str(idx_wire) + ":@")
					file_com.write(atoms_in[idx_bnd])
					file_com.write("\n")

file_com.close()

print("\033[1;32mInstructions printed to '%s'\033[0m" % path_com)
