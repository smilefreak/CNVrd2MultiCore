#!/bin/bash
#SBATCH -J MPI_JOB
#SBATCH -A uoo00009
#SBATCH --time=23:59:00     # Walltime
#SBATCH --ntasks=201        # number of tasks
#SBATCH --mem-per-cpu=4096 # memory/cpu (in MB)

module load R/3.0.3-goolf-1.5.14
mpirun -np 1 R CMD BATCH permutation_mpi.R

