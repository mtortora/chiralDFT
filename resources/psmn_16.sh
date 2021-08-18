#!/bin/bash

#$ -S /bin/bash
#$ -N chiralDFT
#$ -l h_rt=168:00:00
#$ -q CLG6242deb384B,CLG6242deb384C,CLG6226Rdeb192C,CLG6226Rdeb192D,CLG5218deb192C,CLG5218deb192D,CLG5218deb192Th,SLG6142deb384B,SLG6142deb384C
#$ -pe mpi16_debian 256
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
