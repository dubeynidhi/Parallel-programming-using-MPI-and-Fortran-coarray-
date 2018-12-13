#!/bin/bash
#SBATCH -n 1
#SBATCH -J jacobimpi
#SBATCH -o jacobim%j.out
#SBATCH -e jacobim%j.err
module load intel-mpi/2018x
make
time mpirun -n 16 ./jacobimpi