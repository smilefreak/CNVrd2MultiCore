#!/bin/bash

#
# Runs permutations using CNVrd2. 
#
# $1 input_fil

#

cores=10

#@ shell = /bin/bash
#@ environment = COPY_ALL
#@ job_name = mpi_permutation_cnvrd2
#@ job_type = MPICH
#@ total_tasks = 10
#@ account_no = nesi00082
#@ class = merit
#@ group = pd 
#@ blocking = unlimited
#@ resources = ConsumableMemory(2gb) ConsumableVirtualMemory(2gb)
#@ wall_clock_limit = 00:60:00
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

ulimit -v 2097152 -m 2097152
module load R/3.0.1
pwd
Rscript permutation_mpi.R ${cores}

