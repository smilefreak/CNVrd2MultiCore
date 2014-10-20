#!/bin/bash

#
# Runs permutations using CNVrd2. 
#
# $1 input_fil


#@ shell = /bin/bash
#@ environment = COPY_ALL
#@ job_name = mpi_permutation_cnvrd2
#@ job_type = MPICH
#@ total_tasks =201
#@ class = default
#@ account_no = uoo00009
#@ group = nesi
#@ blocking = unlimited
#@ resources = ConsumableMemory(4gb) ConsumableVirtualMemory(4gb)
#@ wall_clock_limit = 23:59:00
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

ulimit -v 4194304 -m 4194304
module load R/3.0.1
cat "$LOADL_HOSTFILE" > hostfile
/usr/mpi/gcc/openmpi-1.4.3/bin/mpirun -hostfile "$LOADL_HOSTFILE" -n 1 /share/apps/R/noarch/3.0.1/gcc-4.4.6/openmpi/1.4.3/bin/R CMD BATCH permutation_mpi.R
