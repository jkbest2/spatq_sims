#!/bin/bash

## Job name
#SBATCH --job-name=fit_densdepq

## Allocations definition
#SBATCH --account=stf
#SBATCH --partition=stf

## Nodes
#SBATCH --nodes=1

## Tasks
#SBATCH --tasks=1

## CPUs per task
#SBATCH --cpus-per-task=28

## Walltime
#SBATCH --time=24:00:00

## Memory per node
#SBATCH --mem=120G

## Set the working directory for this job.
## Directory must exist before submitting job.
#SBATCH --chdir=/gscratch/stf/jkbest/spatq_sims

## Use email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jkbest@uw.edu

export SPATQ_PAR_FITS=7
Rscript 30_fitsims.R densdepq 1 100
