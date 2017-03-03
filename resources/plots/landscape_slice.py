#!/usr/bin/env python3

import sys
import os.path
import argparse
import numpy as np
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Handle command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("path_landscape", help='data grid file')
parser.add_argument("--cmap",         help='matplotlib colormap', default="jet")
parser.add_argument("--x_smooth",     help='x smoothing parameter', type=int, default=1)
parser.add_argument("--y_smooth",     help='y smoothing parameter', type=int, default=1)
parser.add_argument("--lines",        help='display iso-y lines', dest="lines", action="store_true")
parser.add_argument("--normalize",    help='normalise landscape', dest="normalize", action="store_true")

args = parser.parse_args()

parser.set_defaults(lines=False)

# Load data
if not os.path.isfile(args.path_landscape):
	print("\033[1;31mCouldn't read landscape file %s - aborting\033[0m" % args.path_landscape)
	sys.exit()

with open(args.path_landscape, mode="r") as file_landscape:
	for idx,line in enumerate(file_landscape):
		if not line.rstrip(): break
	n_x  = idx

data_landscape = np.genfromtxt(args.path_landscape)

n_y            = len(data_landscape) // n_x

x_smooth       = args.x_smooth
y_smooth       = args.y_smooth

if n_x % x_smooth != 0:
	print("\033[1;31mx_smooth must divide %d - choose a different value\033[0m" % n_x)
	sys.exit()

if n_y % y_smooth != 0:
	print("\033[1;31my_smooth must divide %d - choose a different value\033[0m" % n_y)
	sys.exit()

data_landscape              = data_landscape.reshape(n_y, n_x, 3)

n_x_smooth                  = n_x // x_smooth
n_y_smooth                  = n_y // y_smooth

# Rolling averages over r and theta
data_landscape              = data_landscape.reshape(n_y_smooth, y_smooth, n_x, 3)
data_landscape              = np.mean(data_landscape, axis=1)

data_landscape              = data_landscape.reshape(n_y_smooth, n_x_smooth, x_smooth, 3)
data_landscape              = np.mean(data_landscape, axis=2)

# Setup colormap for r
y_min                       = np.min(data_landscape[..., 0])
y_max                       = np.max(data_landscape[..., 0])

cmap                        = plt.get_cmap(args.cmap)
cNorm                       = colors.Normalize(vmin=y_min, vmax=y_max)

scalarMap                   = cmx.ScalarMappable(cmap=cmap, norm=cNorm)
scalarMap._A                = []

# Switch x-range from [0:pi] to [-pi/2:pi/2]
data                        = np.zeros_like(data_landscape)

if n_x_smooth % 2 == 0:
	data[:, :n_x_smooth//2, :] = data_landscape[:, n_x_smooth//2:, :].copy()
	data[:, n_x_smooth//2:, :] = data_landscape[:, :n_x_smooth//2, :].copy()

else:
	data[:, :n_x_smooth//2, :] = data_landscape[:, n_x_smooth//2+1:, :].copy()
	data[:, n_x_smooth//2:, :] = data_landscape[:, :n_x_smooth//2+1, :].copy()

data[:, :n_x_smooth//2, 1] -= np.pi

if args.normalize: data[...,2] /= np.max(np.abs(data[...,2]), axis=1, keepdims=True)

# Setup plot window
plt.style.use('seaborn-whitegrid')
plt.rc('font', **{'family':'serif','size':'14'})

plt.figure(figsize=(9,6))

plt.xlim([-np.pi/2, np.pi/2])
plt.ylim([np.min(data[...,2]), np.max(data[...,2])])

for i in range(n_y_smooth):
	if args.lines: plt.plot(data[i,:,1], data[i,:,2], color=scalarMap.to_rgba(data[i,0,0]))
	else:          plt.fill_between(data[i,:,1], 0, data[i,:,2], color=scalarMap.to_rgba(data[i,0,0]), alpha=0.25, lw=0)

cb = plt.colorbar(scalarMap)
				  
plt.show()
