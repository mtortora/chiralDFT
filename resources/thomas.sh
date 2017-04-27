#!/bin/bash -l

#$ -S /bin/bash

#$ -l h_rt=40:00:0
#$ -N cDFT-DNA
#$ -pe mpi 864

#$ -wd TARGETPATH/../

module unload compilers/intel/2017/update1
module unload mpi/intel/2017/update1/intel

module load compilers/gnu/4.9.2
module load mpi/openmpi/1.10.1/gnu-4.9.2 

gerun TARGETPATH/TARGET > DATPATH/log.out
