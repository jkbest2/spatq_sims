#!/bin/bash

## Job name
#SBATCH --job-name=prep_sims

## Allocations definition
#SBATCH --account=stf
#SBATCH --partition=stf

## Nodes
#SBATCH --nodes=1

## Walltime (2 hours)
#SBATCH --time=2:00:00

## Number of cores per node
#SBATCH --ntasks-per-node=28

## Memory per node
#SBATCH --mem=120G

## Set the working directory for this job.
## Directory must exist before submitting job.
#SBATCH --chdir=/gscratch/stf/jkbest/spatq_sims

## Use email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jkbest@uw.edu

julia --project -- 10_prepsims.jl
