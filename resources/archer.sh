#!/bin/bash --login

#PBS -N s6hb-n
#PBS -l select=128
#PBS -l walltime=24:00:00
#PBS -m bee
#PBS -M maxime.tortora@icloud.com
#PBS -A e280-Doye


module switch PrgEnv-cray/5.2.56 PrgEnv-gnu/5.2.56
module switch gcc/5.1.0 gcc/6.1.0

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

aprun -n 3072 TARGETPATH/TARGET > DATPATH/log.out
