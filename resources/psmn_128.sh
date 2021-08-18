#!/bin/bash

#$ -S /bin/bash
#$ -N Y21M
#$ -l h_rt=168:00:00
#$ -q Epyc7702deb512
#$ -pe mpi128_debian 640
#$ -cwd
#$ -V
#$ -m be

module load GCC/7.2.0/OpenMPI/3.0.0

# given by SGE
HOSTFILE="${TMPDIR}/machines"

cd "${SGE_O_WORKDIR}" || { echo "cannot cd to ${SGE_O_WORKDIR}"; exit 1; }

PREFIX="/applis/PSMN/debian9/software/Compiler/GCC/7.2.0/OpenMPI/3.0.0"
MPIRUN="${PREFIX}/bin/mpirun"

"${MPIRUN}" -prefix "${PREFIX}" -hostfile "${HOSTFILE}" -np "${NSLOTS}" TARGETPATH/TARGET > DATPATH/log.out
