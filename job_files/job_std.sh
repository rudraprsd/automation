#/bin/sh
#PBS -N Si_scf
#PBS -o mpi-out
#PBS -e mpi-error
#PBS -l select=1:ncpus=16
#PBS -q workq 
cd $PBS_O_WORKDIR


export PATH=/opt/apps/gcc/gcc-11.2.0/bin:$PATH
LD_LIBRARY_PATH=/opt/apps/gcc/gcc-11.2.0/lib64:$PATH
source /opt/apps/oneapi/2023.2.0/setvars.sh --force
export PATH=/opt/apps/vasp/vasp610/w90v12/bin:$PATH

mpirun -np 16 vasp_std >out

