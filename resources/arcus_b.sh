#!/bin/bash

#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --time=120:00:00
#SBATCH --partition=compute
#SBATCH --job-name=cDFT-DNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxime.tortora@icloud.com


module load mvapich2/2.1.0__gcc-4.9.2

. enable_arcus-b_mpi.sh

mpirun $MPI_HOSTS TARGETPATH/TARGET &> DATPATH/$SLURM_JOB_ID.out
