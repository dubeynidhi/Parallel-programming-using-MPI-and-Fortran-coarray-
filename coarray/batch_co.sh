#!/bin/bash
#SBATCH -n 1
#SBATCH -J jacobicoarray
#SBATCH -o jacobi%j.out
#SBATCH -e jacobi%j.err
module load intel-mpi/2018x
make
export FOR_COARRAY_NUM_IMAGES=16
time ./jacobicoarray 