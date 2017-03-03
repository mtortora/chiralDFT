#!/usr/bin/env python3

import sys
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from cycler import cycler


# Handle command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("path_landscape", help='data grid file')
parser.add_argument("--cmap",         help='matplotlib colormap', default="jet")
parser.add_argument("--normalize",    help='normalise landscape', dest="normalize", action="store_true")
parser.add_argument("--show_min",     help='display minima',      dest="show_min",  action="store_true")

parser.set_defaults(normalize=False)
parser.set_defaults(show_min=False)

args = parser.parse_args()

# Load data
if not os.path.isfile(args.path_landscape):
	print("\033[1;31mCouldn't read landscape file %s - aborting\033[0m" % args.path_landscape)
	sys.exit()

with open(args.path_landscape, mode="r") as file_landscape:
	for idx,line in enumerate(file_landscape):
		if not line.rstrip(): break
	ny = idx

data_landscape = np.genfromtxt(args.path_landscape)

nx             = len(data_landscape) // ny

X              = data_landscape[::ny,0]
Y              = data_landscape[:ny,1]
Z              = data_landscape[:,2]

T              = Z.reshape(nx,ny).T

if args.normalize: T /= np.max(np.abs(T), axis=0, keepdims=True)

# Colorbar font formatter
def fmt(x, pos):
	a, b = '{:.0e}'.format(x).split('e')
	b = int(b)
	if x != 0: return r'${} \times 10^{{{}}}$'.format(a, b)
	else: return r'$0$'

# Setup plot window
fontdict   = {'family':'serif','size':'14'}
facecolors = ["#FBC15E","#8EBA42","#988ED5","#348ABD","#FFB5B8","#E24A33","#777777"]

plt.style.use('ggplot')
plt.rc('font',**fontdict)
plt.rcParams['axes.prop_cycle'] = cycler(color=facecolors)

for key in ['fancybox','shadow']: plt.rcParams['legend.'+key] = True

plt.figure(figsize=(10,6))

# Plot landscape and minima
if args.show_min:
	data_min      = np.zeros([nx, 2])
	idx_min       = np.argmin(T, axis=0)
	
	data_min[:,0] = X
	data_min[:,1] = Y[idx_min]
	
	plt.plot(data_min[:,0], data_min[:,1], "k--", lw=3)

implot = plt.imshow(T, aspect='auto', origin='lower', extent=(X.min(),X.max(),Y.min(),Y.max()))
implot.set_cmap(args.cmap)

plt.colorbar(format=ticker.FuncFormatter(fmt))

plt.show()
