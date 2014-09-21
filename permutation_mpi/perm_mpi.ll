#!/bin/bash

#
# Runs permutations using CNVrd2. 
#
# $1 input_fil

#

cores=200

#@ shell = /bin/bash
#@ environment = COPY_ALL
#@ job_name = mpi_permutation_cnvrd2
#@ job_type = MPICH
#@ total_tasks = 200
#@ class = default
#@ group = nesi
#@ blocking = unlimited
#@ resources = ConsumableMemory(2gb) ConsumableVirtualMemory(2gb)
#@ wall_clock_limit = 23:59:00
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

ulimit -v 2097152 -m 2097152
module load R/3.0.1
Rscript permutation_mpi.R ${cores}

