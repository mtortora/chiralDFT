#!/usr/bin/env python

import numpy as np
import os.path
import struct
import sys


# Input manager
if len(sys.argv) != 3:
	print("\033[1;31mUsage is %s data_folder n_reduce\033[0m" % sys.argv[0])
	sys.exit()

path_data = sys.argv[1].rstrip("/")
n_reduce  = int(sys.argv[2])

if os.path.isdir(path_data):
	path_ref = path_data + '/q_pert.out'
	path_kt  = path_data + '/kt_threaded.out'
	path_k2  = path_data + '/k2_threaded.out'

	for path in [path_ref, path_kt, path_k2]:
		if not os.path.isfile(path):
			print("\033[1;31mCouldn't read file %s - aborting\033[0m" % path)
			sys.exit()

else:
	print("\033[1;31mCouldn't open directory %s - aborting\033[0m" % path_data)
	sys.exit()


# Read datafiles
data_ref    = np.genfromtxt(path_ref)
eta_grid    = data_ref[:, 0]
n_steps_eta = len(eta_grid)

data        = []

for path in [path_kt, path_k2]:
	with open(path, mode="rb") as file: file_content = file.read()

	file_data = struct.unpack("d" * (len(file_content) // 8), file_content)
	data.append(file_data)

# Format data
data          = np.asarray(data)
n_threads_in  = np.size(data) // (2*n_steps_eta)
n_threads_out = n_threads_in  // n_reduce

if n_threads_in % n_reduce != 0:
	print("\033[1;31mn_reduce needs to divide n_threads (%d)\033[0m" % n_threads_in)
	sys.exit()

# Process local averages
data   = np.reshape(data, (2, n_steps_eta, n_threads_out, n_reduce), order="F")
kt, k2 = np.mean(data, axis=-1)
q      = kt / k2

# Setup output containers
kt_res = np.zeros_like(data_ref)
k2_res = np.zeros_like(data_ref)
q_res  = np.zeros_like(data_ref)

kt_res[:, 0] = eta_grid
kt_res[:, 1] = np.mean(kt, axis=-1)
kt_res[:, 2] = np.min (kt, axis=-1)
kt_res[:, 3] = np.max (kt, axis=-1)

k2_res[:, 0] = eta_grid
k2_res[:, 1] = np.mean(k2, axis=-1)
k2_res[:, 2] = np.min (k2, axis=-1)
k2_res[:, 3] = np.max (k2, axis=-1)

q_res[:, 0]  = eta_grid
q_res[:, 1]  = np.mean(q,  axis=-1)
q_res[:, 2]  = np.min (q,  axis=-1)
q_res[:, 3]  = np.max (q,  axis=-1)

# Save to files
np.savetxt(path_data + "/kt.res", kt_res)
np.savetxt(path_data + "/k2.res", k2_res)
np.savetxt(path_data + "/q_pert.res", q_res)