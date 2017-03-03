#!/usr/bin/env python

import os
import sys
import subprocess

import numpy as np


# Input/output
if len(sys.argv) != 5:
	print("\033[1;31mUsage is %s trajectory topology idx_conf num_box\033[0m" % sys.argv[0])
	sys.exit()


POS_MM_BACK1 = -0.3400
POS_MM_BACK2 = 0.3408
POS_BACK     = -0.4
POS_BASE     = 0.4

padding      = 0.
radius       = 0.1

color        = 'r="0.2734" g="0.5078" b="0.7031"'


def PCA(points_in):
	# Translate center-of-mass back to the origin
	points_out  = np.asarray(points_in)
	cm          = np.mean(points_out, axis=0)
	
	points_out -= cm
	points_out  = points_out.T
	
	# Perform Principal Component Analysis of the configuration by singular-value decomposition
	P, D, Q     = np.linalg.svd(points_out)
	
	# P is the rotation matrix expressing the covariance matrix eigenvectors in the reference frame
	rot         = np.roll(P, 2, axis=1)
	inv_rot     = rot.T
	
	points_out  = np.dot(inv_rot, points_out)
	
	return cm, rot, points_out


path_conf    = os.path.splitext(sys.argv[1])[0] + ".cnf"
path_box     = os.path.splitext(sys.argv[1])[0] + ".cmm"

path_traj    = sys.argv[1]
path_top     = sys.argv[2]


file_traj    = open(path_traj, mode="r")
file_top     = open(path_top,  mode="r")
file_conf    = open(path_conf, mode="w")
file_output  = open(path_box,  mode="w")

idx_conf     = int(sys.argv[3])
num_box      = int(sys.argv[4])

conf_counter = 0

file_top.readline()

nucleotides  = file_top.readlines()
timeline     = file_traj.readline()

while timeline:
	conf_counter += 1
	
	line_box      = file_traj.readline()
	line_E        = file_traj.readline()
	
	if conf_counter == idx_conf:
		file_conf.write(timeline)
		file_conf.write(line_box)
		file_conf.write(line_E)
		
		cm_pos = []

		for nucleotide in nucleotides:
			nucleo_line = file_traj.readline()
			nucleo_data = nucleo_line.split()
		
			ci   = [float(x) for x in nucleo_data[0:3]]
			a1   = [float(x) for x in nucleo_data[3:6]]
			a3   = [float(x) for x in nucleo_data[6:9]]
	
			ci   = np.asarray(ci)
			a1   = np.asarray(a1)
			a3   = np.asarray(a3)
			a2   = np.cross(a3,a1)
		
			base = ci + a1*POS_BASE
			back = ci + a1*POS_MM_BACK1 + a2*POS_MM_BACK2
			#back = ci + a1*POS_BACK
		
			file_conf.write(nucleo_line)

			cm_pos.append(back)

		break

	else:
		for nucleotide in nucleotides: file_traj.readline()

	timeline = file_traj.readline()


if conf_counter < idx_conf:
	print("\033[1;31mCould only find %d configuration(s) in trajectory file\033[0m" % conf_counter)
	sys.exit()


cm, rot, points = PCA(cm_pos)

boxes           = np.zeros([num_box, 3, 8])
subdiv          = np.linspace(np.min(points[2, :]), np.max(points[2, :]), num=num_box+1)

file_output.write('<marker_set name="marker set 1">\n')

for idx_box, box in enumerate(boxes):
	vtx_box = []

	for point in points.T:
		if subdiv[idx_box] <= point[2] <= subdiv[idx_box+1]: vtx_box.append(point)

	cm_box, rot_box, points_box = PCA(vtx_box)
	
	center    = (np.max(points_box, axis=1) + np.min(points_box, axis=1))/2.
	size      = (np.max(points_box, axis=1) - np.min(points_box, axis=1))/2. + padding

	box[0, 0] = -size[0]
	box[1, 0] = -size[1]
	box[2, 0] = -size[2]

	box[0, 1] = +size[0]
	box[1, 1] = -size[1]
	box[2, 1] = -size[2]

	box[0, 2] = +size[0]
	box[1, 2] = +size[1]
	box[2, 2] = -size[2]

	box[0, 3] = -size[0]
	box[1, 3] = +size[1]
	box[2, 3] = -size[2]

	box[0, 4] = -size[0]
	box[1, 4] = -size[1]
	box[2, 4] = +size[2]

	box[0, 5] = +size[0]
	box[1, 5] = -size[1]
	box[2, 5] = +size[2]

	box[0, 6] = +size[0]
	box[1, 6] = +size[1]
	box[2, 6] = +size[2]

	box[0, 7] = -size[0]
	box[1, 7] = +size[1]
	box[2, 7] = +size[2]

	box      += center[:, np.newaxis]
	
	box       = np.dot(rot_box, box) + cm_box[:, np.newaxis]
	box       = np.dot(rot,     box) #+ cm    [:, np.newaxis]

	links     = [[0,1], [1,2], [2,3], [3,0], [4,5], [5,6], [6,7], [7,4], [0,4], [1,5], [2,6], [3,7]]

	for idx_marker, marker in enumerate(box.T):
		idx = 8*idx_box + idx_marker
		file_output.write('<marker id ="' + str(idx+1) + '" ')
		file_output.write('x="' + str(marker[0]) + '" ')
		file_output.write('y="' + str(marker[1]) + '" ')
		file_output.write('z="' + str(marker[2]) + '" ')
		file_output.write(color + ' ')
		file_output.write('radius="' + str(radius) + '"/>')
		file_output.write('\n')

	for link in links:
		file_output.write('<link ')
		file_output.write('id1="' + str(8*idx_box+1 + link[0]) + '" ')
		file_output.write('id2="' + str(8*idx_box+1 + link[1]) + '" ')
		file_output.write(color + ' ')
		file_output.write('radius="' + str(radius) + '"/>')
		file_output.write('\n')

file_output.write('</marker_set>')

file_traj.close()
file_top.close()
file_conf.close()
file_output.close()

print("\033[1;32mConfiguration printed to '%s'\033[0m" % path_conf)
print("\033[1;32mBox vertices printed to '%s'\033[0m" % path_box)

path = os.environ.get("PATH_TO_OXDNA")

if path:
	os.system(path + "/UTILS/traj2chimera.py " + path_conf + " " + path_top)
	
	sp = subprocess.Popen(["/bin/bash", "-i", "-c", "chimera " + path_conf + ".pdb" + " chimera.com " + path_box])
	sp.communicate()
