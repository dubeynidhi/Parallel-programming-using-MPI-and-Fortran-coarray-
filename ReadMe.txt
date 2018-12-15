This is a Jacobian iteration implementation on a 4*4 processor arrangement using MPI and Fortran coarray. Each process is implemented as a 66*66 matrix and true solution is implemented as 256*256. The start and end index of inner elements for both u matrix and f and t matrix process wise. The update the values, transfer data and check for convergence. Finally, the covergence value and the final difference is calculated and printed on screen in root process or image 1. For detailed in description of the project refer the problem formulation pdf.


The processors are in this form for MPI

0	1	2	3
4	5	6	7
8	9	10	11
12	13	14	15

The processors are in this form for coarray
1	2	3	4
5	6	7	8
9	10	11	12
13	14	15	16



how to run:
Load the modules for Fortran and MPI and then:

For MPI:
type: sbatch batch_mpi.sh
You will see a message with your job number on the screen
.out file will create the output and .err will create the errors
To clean, type make clean

For Coarray:
type: sbatch batch_co.sh
You will see a message with your job number on the screen
.out file will create the output and .err will create the errors
To clean, type make clean

