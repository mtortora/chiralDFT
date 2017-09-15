#!/bin/bash --login

#PBS -N s-6ho
#PBS -l select=128:bigmem=true
#PBS -l walltime=24:00:00
#PBS -m bee
#PBS -M maxime.tortora@icloud.com
#PBS -A e552


module switch PrgEnv-cray/5.2.82 PrgEnv-gnu/5.2.82

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

aprun -n 3072 TARGETPATH/TARGET > DATPATH/log.out
