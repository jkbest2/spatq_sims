#!/bin/bash

## Job name
#SBATCH --job-name=run_sims

## Allocations definition
#SBATCH --account=stf
#SBATCH --partition=stf

## Nodes
#SBATCH --nodes=1

## Cores per node
#SBATCH --ntasks-per-node=28

## Walltime (2 hours)
#SBATCH --time=1:00:00

## Memory per node
#SBATCH --mem=100G

## Set the working directory for this job.
## Directory must exist before submitting job.
#SBATCH --chdir=/gscratch/stf/jkbest/spatq_sims

## Use email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jkbest@uw.edu

module load contrib/julia/1.3

julia --project -p 28 --check-bounds=no -L 20_gensims.jl 21_runsims.jl

