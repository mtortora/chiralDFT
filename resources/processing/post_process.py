#!/usr/bin/env python3

import numpy as np
import os.path
import struct
import sys

from numpy.polynomial.legendre import Legendre


if sys.version_info < (3, 0):
	print("\033[1;31mPython 2 or older detected - switch to Python 3\033[0m")
	sys.exit()

if len(sys.argv) != 6:
	print("\033[1;31mUsage is %s data_folder file_globals eta_min eta_max n_reduce\033[0m" % sys.argv[0])
	sys.exit()


# ODF parameters
n_steps_eta   = 200
n_steps_theta = 200
n_steps_max   = 1000

tol_odf       = 1e-6

path_data     = sys.argv[1].rstrip("/")
path_globals  = sys.argv[2]

eta_min       = float(sys.argv[3])
eta_max       = float(sys.argv[4])

n_reduce      = int(sys.argv[5])

parameters    = ["Q_MIN","Q_MAX","N_STEPS_Q","N_L","N_STEPS_THETA","MODE"]
param_dict    = {}

# Read parameters from globals.hpp
if os.path.isfile(path_globals):
	with open(path_globals, mode="r") as file_globals:
		for line in file_globals:
			line_data = line.split()
			
			if len(line_data) == 3:
				param = line_data[1]
				if param in parameters: param_dict[param] = line_data[2]

	for param in parameters[:2]: param_dict[param] = float(param_dict[param])
	for param in parameters[2:]: param_dict[param] = int  (param_dict[param])

else:
	print("\033[1;31mCouldn't read file %s - aborting\033[0m" % path_globals)
	sys.exit()

# Read simulation data
if os.path.isdir(path_data):
	path_ref = path_data + '/excluded_matrix.out'
	
	path_k2  = path_data + '/k2_threaded.out'
	path_kt  = path_data + '/kt_threaded.out'
	
	eta_grid = np.genfromtxt(path_data + '/k2.out')[:,0]
	n_eta    = len(eta_grid)

	K        = []
	
	# Local average of threaded K2 and kt
	for path in [path_k2, path_kt]:
		if os.path.isfile(path):
			with open(path, mode="rb") as file: file_content = file.read()

			data          = struct.unpack("d" * (len(file_content) // 8), file_content)
			data          = np.asarray(data)
		
			n_threads_in  = np.size(data) // n_eta
			n_threads_out = n_threads_in  // n_reduce

			data_res      = np.zeros([n_eta, n_threads_out])

			# Partial reduction over n_reduce threads
			for idx_thread_out in range(n_threads_out):
				data_sum = np.zeros(n_eta)
		
				for idx_reduce in range(n_reduce):
					idx_thread_in = idx_thread_out*n_reduce + idx_reduce
				
					data_loc      = data[idx_thread_in * n_eta:(idx_thread_in+1) * n_eta]
					data_sum     += data_loc / n_reduce
			
				data_res[:, idx_thread_out] = data_sum
		
		else:
			print("\033[1;31mCouldn't read file %s - aborting\033[0m" % path)
			sys.exit()

		K.append(data_res)
			 
	k2, kt      = K
	q           = kt/k2
	
	q_ave       = np.mean(q, axis=1)
	q_min       = np.min (q, axis=1)
	q_max       = np.max (q, axis=1)
	
	data_q_pert = np.zeros([n_eta, 4])
	
	data_q_pert[:,0] = eta_grid
	data_q_pert[:,1] = q_ave
	data_q_pert[:,2] = q_min
	data_q_pert[:,3] = q_max

	path_q_pert      = path_data + '/q_pert.res'
	np.savetxt(path_q_pert, data_q_pert)

	print("\033[1;32mPerturbative pitch printed to '%s'\033[0m" % path_q_pert)

	# Read particle volume from legendre_matrix.out
	if os.path.isfile(path_ref):
		with open(path_ref, mode="r") as file_ref:
			V0   = float(file_ref.readline().split()[2])
			Veff = float(file_ref.readline().split()[2])
		excluded_ref = np.genfromtxt(path_ref, skip_header=2)

		if   len(excluded_ref) == param_dict["N_STEPS_THETA"]: mode = "odf_full"
		elif len(excluded_ref) == param_dict["N_L"]:           mode = "odf_legendre"
		else:
			print("\033[1;31m%s is incompatible with ODF parameters - aborting\033[0m" % path_ref)
			sys.exit()
	
	else:
		print("\033[1;31mCouldn't read file %s - aborting\033[0m" % path_ref)
		sys.exit()

else:
	print("\033[1;31mCouldn't open directory %s - aborting\033[0m" % path_data)
	sys.exit()


# Binary file converter
def LoadBinary(Q_MIN, Q_MAX, N_STEPS_Q, N_L, **kwargs):
	q_grid        = np.linspace(Q_MIN, Q_MAX, N_STEPS_Q)
	data_res      = []
	
	for q in q_grid:
		q_form        = "{:.6f}".format(q)
		path_legendre = path_data + '/legendre_matrix_' + q_form + '.out'
		
		if not os.path.isfile(path_legendre):
			print("\033[1;31mCouldn't read file %s - check range parameters\033[0m" % path_legendre)
			sys.exit()
		
		with open(path_legendre, mode="rb") as file_legendre: file_content = file_legendre.read()

		data          = struct.unpack("d" * (len(file_content) // 8), file_content)
		data          = np.asarray(data)

		n_threads_in  = np.size(data) // N_L**2
		n_threads_out = n_threads_in  // n_reduce

		# Partial reduction over n_reduce threads
		for idx_thread_out in range(n_threads_out):
			data_sum = np.zeros([N_L, N_L])
			
			for idx_reduce in range(n_reduce):
				idx_thread_in = idx_thread_out*n_reduce + idx_reduce
				
				data_loc      = data[idx_thread_in * N_L**2:(idx_thread_in+1) * N_L**2]
				data_loc      = data_loc.reshape(N_L, N_L, order="F")
				data_sum     += data_loc / n_reduce
			
			data_res.append(data_sum)

	data_res = np.asarray(data_res)
	data_res = data_res.reshape(N_STEPS_Q, n_threads_out, N_L, N_L, order="C")

	return data_res


# ODF optimiser from angle-dependant excluded volume
def ODFOnsagerOptimiser(E):
	converged  = False
	d_theta    = np.pi / (param_dict["N_STEPS_THETA"]-1.)
	
	eta_grid   = np.linspace(eta_min, eta_max, n_steps_eta)
	theta_grid = np.linspace(0., np.pi, param_dict["N_STEPS_THETA"])
	
	S_res      = np.zeros(n_steps_eta)
	P_res      = np.zeros(n_steps_eta)
	F_tot      = np.zeros(n_steps_eta)

	path_p     = path_data + '/p.res'
	path_psi   = path_data + '/psi.res'
	path_s     = path_data + '/order_param.res'
	
	with  open(path_psi, mode="wb") as file_psi:
		for idx_eta, eta in enumerate(eta_grid):
			n_dens  = eta / V0
			eta_eff = eta * Veff/V0
			g_pl    = (1-0.75*eta_eff) / (1-eta_eff)**2
			
			psi     = np.exp(-theta_grid**2)
			
			# ODF self-consistent solver
			for idx_iter in range(n_steps_max):
				psi_dummy = psi.copy()
				
				psi       = np.exp(-g_pl * n_dens/(2*np.pi)**2 * ((E + E.T)/2. * np.sin(theta_grid) * psi_dummy).sum(axis=1) * d_theta)
				psi      /= (2*np.pi)**2 * d_theta * (psi * np.sin(theta_grid)).sum()
				
				if np.max(np.abs(psi - psi_dummy)) < tol_odf:
					converged = True
					break
		
			if not converged:
				print("\033[1;31mConvergence failed - eta=%f\033[0m" % eta)
				sys.exit()
	
			# Free energy
			psi_p = psi*np.sin(theta_grid)
			b2    = (np.outer(psi_p, psi_p) * E).sum() * d_theta**2/2.
					
			f_id  = n_dens * (2*np.pi)**2 * d_theta * (np.sin(theta_grid) * psi * np.log(psi)).sum()
			f_exc = b2/V0 * n_dens * eta*g_pl
			
			p_id  = eta
			p_exc = b2/V0 * eta * (eta_eff - eta_eff**2/2.) / (1-eta_eff)**3
			
			F_tot[idx_eta] = f_id + f_exc
			P_res[idx_eta] = p_id + p_exc
			
			# Order parameter
			s_tmp          = (3*np.cos(theta_grid)**2 - 1) / 2.
			s_tmp         *= (2*np.pi)**2 * d_theta * np.sin(theta_grid) * psi
			
			S_res[idx_eta] = s_tmp.sum()

			# Save ODF
			data_odf       = np.zeros([len(theta_grid), 3])
	
			data_odf[:, 0] = eta * np.ones(len(theta_grid))
			data_odf[:, 1] = theta_grid
			data_odf[:, 2] = psi
				
			np.savetxt(file_psi, data_odf)
			file_psi.write(bytes('\n', "UTF-8"))

	data_p       = np.zeros([n_steps_eta, 2])
	data_p[:, 0] = eta_grid
	data_p[:, 1] = P_res
	
	data_s       = np.zeros([n_steps_eta, 2])
	data_s[:, 0] = eta_grid
	data_s[:, 1] = S_res

	np.savetxt(path_p, data_p)
	np.savetxt(path_s, data_s)

	print("\033[1;32mODFs printed to '%s'\033[0m" % path_psi)
	print("\033[1;32mOrder parameter printed to '%s'\033[0m" % path_s)
	print("\033[1;32mOsmotic pressure printed to '%s'\033[0m" % path_p)
	
	return F_tot,


# ODF optimiser from Legendre-projected excluded volume
def ODFLegendreOptimiser(data_q, reference_run):
	n_threads     = np.shape(data_q)[0]
	n_l           = np.shape(data_q)[-1]
	
	d_theta       = np.pi / (n_steps_theta-1.)
	
	eta_grid      = np.linspace(eta_min, eta_max, n_steps_eta)
	theta_grid    = np.linspace(0., np.pi, n_steps_theta)
	
	coeffs        = np.sqrt(np.arange(n_l) + 0.5)
	
	S_res         = np.zeros(n_steps_eta)
	F_tot         = np.zeros([n_threads, n_steps_eta])
	
	for idx_thread,E_q_loc in enumerate(data_q):
		print("Processing thread %d out of %d" % tuple([idx_thread+1, n_threads]))
		sys.stdout.flush()
		
		converged = False
		
		for idx_eta, eta in enumerate(eta_grid):
			n_dens = eta / V0
			g_pl   = (1-0.75*eta) / (1-eta)**2
			
			psi    = np.exp(-theta_grid**2)
			psi_l  = np.zeros(n_l)

			# ODF self-consistent solver
			for idx_iter in range(n_steps_max):
				for l in range(n_l):
					if np.any(E_q_loc[l,:] != 0.):
						order    = np.zeros(l+1)
						order[l] = 1
						psi_l[l] = (Legendre(order)(np.cos(theta_grid)) * psi * np.sin(theta_grid)).sum()

				psi_l      *= coeffs * d_theta
				psi_dummy   = psi.copy()

				coeffs_excl = np.dot(E_q_loc+E_q_loc.T, psi_l) / 2.

				psi         = np.exp(-g_pl * n_dens/(2*np.pi)**2 * Legendre(coeffs * coeffs_excl)(np.cos(theta_grid)))
				psi        /= (2*np.pi)**2 * d_theta * (psi * np.sin(theta_grid)).sum()
									 
				if np.max(np.abs(psi - psi_dummy)) < tol_odf:
					converged = True
					break

			if not converged:
				print("\033[1;31mConvergence failed - eta=%f\033[0m" % eta)
				sys.exit()

			# Free energy
			b2    = (np.outer(psi_l, psi_l) * E_q_loc).sum() / 2.

			f_id  = n_dens * (2*np.pi)**2 * d_theta * (np.sin(theta_grid) * psi * np.log(psi)).sum()
			f_exc = b2/V0 * n_dens * eta*g_pl

			F_tot[idx_thread, idx_eta] = f_id + f_exc

			# Order parameter
			if reference_run:
				s_tmp          = (3*np.cos(theta_grid)**2 - 1) / 2.
				s_tmp         *= (2*np.pi)**2 * d_theta * np.sin(theta_grid) * psi

				S_res[idx_eta] = s_tmp.sum()

		if reference_run:
			data_res       = np.zeros([n_steps_eta, 2])
			data_res[:, 0] = eta_grid
			data_res[:, 1] = S_res

			path_order     = path_data + '/order_param.res'
			np.savetxt(path_order, data_res)

			print("\033[1;32mOrder parameter printed to '%s'\033[0m" % path_order)
			
	return F_tot


# Preliminary run
print("\033[1;34mAggregated reference run\033[0m")

if mode == "odf_legendre":
	data_ref = excluded_ref[np.newaxis,...]
	F_ref    = ODFLegendreOptimiser(data_ref, True)

if mode == "odf_full": F_ref = np.asarray(ODFOnsagerOptimiser(excluded_ref))

# Chiral landscape sweeping run
if param_dict["MODE"] == 0:
	data      = LoadBinary(**param_dict)

	q_min     = param_dict["Q_MIN"]
	q_max     = param_dict["Q_MAX"]

	n_steps_q = np.shape(data)[0]
	n_threads = np.shape(data)[1]

	eta_grid  = np.linspace(eta_min, eta_max, n_steps_eta)
	q_grid    = np.linspace(q_min,   q_max,   n_steps_q)

	F_tot     = np.zeros([n_steps_q, n_threads, n_steps_eta])
	
	for idx_q,data_q in enumerate(data):
		print("\033[1;34mMacroscopic pitch %d out of %d\033[0m" % tuple([idx_q+1, n_steps_q]))
		
		F_tot[idx_q,...] = ODFLegendreOptimiser(data_q, False)

	F_res              = np.mean(F_tot, axis=1)

	idx_min            = np.argmin(F_tot, axis=0)
	q_eq               = q_grid[idx_min]

	# Work out K constants as successive differentials of the free energy
	dq                 = np.diff(q_grid)[0]

	Kt                 = np.diff(F_tot, axis=0, n=1) / dq
	K2                 = np.diff(F_tot, axis=0, n=2) / dq**2

	# Average Kt, K2 over all q's
	Kt                 = np.mean(Kt, axis=0)
	K2                 = np.mean(K2, axis=0)

	q_res              = np.mean(q_eq, axis=0)
	q_inf              = np.min (q_eq, axis=0)
	q_sup              = np.max (q_eq, axis=0)

	Kt_res             = np.mean(Kt, axis=0)
	Kt_inf             = np.min (Kt, axis=0)
	Kt_sup             = np.max (Kt, axis=0)

	K2_res             = np.mean(K2, axis=0)
	K2_inf             = np.min (K2, axis=0)
	K2_sup             = np.max (K2, axis=0)

	# Format data for 3d visualisation
	X, Y               = np.meshgrid(eta_grid, q_grid)
	Z                  = F_res - F_ref

	data_f_res         = np.zeros([n_steps_eta, n_steps_q, 3])
	data_q_res         = np.zeros([n_steps_eta, 4])

	data_kt_res        = np.zeros([n_steps_eta, 4])
	data_k2_res        = np.zeros([n_steps_eta, 4])

	data_f_res[..., 0] = X.T
	data_f_res[..., 1] = Y.T
	data_f_res[..., 2] = Z.T

	data_q_res[:, 0]   = eta_grid
	data_q_res[:, 1]   = q_res
	data_q_res[:, 2]   = q_inf
	data_q_res[:, 3]   = q_sup

	data_kt_res[:, 0]  = eta_grid
	data_kt_res[:, 1]  = Kt_res
	data_kt_res[:, 2]  = Kt_inf
	data_kt_res[:, 3]  = Kt_sup

	data_k2_res[:, 0]  = eta_grid
	data_k2_res[:, 1]  = K2_res
	data_k2_res[:, 2]  = K2_inf
	data_k2_res[:, 3]  = K2_sup

	path_landscape     = path_data + '/energy_landscape.res'
	path_q             = path_data + '/qmin.res'

	path_kt            = path_data + '/kt.res'
	path_k2            = path_data + '/k2.res'

	with open(path_landscape, mode="wb") as file_landscape:
		for data_slice in data_f_res:
			np.savetxt(file_landscape, data_slice)
			file_landscape.write(bytes('\n', "UTF-8"))

	np.savetxt(path_q, data_q_res)

	np.savetxt(path_kt, data_kt_res)
	np.savetxt(path_k2, data_k2_res)

	print("\033[1;32mLandscape printed to '%s'\033[0m" % path_landscape)
	print("\033[1;32mPitch printed to '%s'\033[0m" % path_q)
	print("\033[1;32mKt printed to '%s'\033[0m" % path_kt)
	print("\033[1;32mK2 printed to '%s'\033[0m" % path_k2)
