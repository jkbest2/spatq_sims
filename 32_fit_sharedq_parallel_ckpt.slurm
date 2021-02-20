#!/bin/bash

## Job name
#SBATCH --job-name=fit_sharedq

## Allocations definition
#SBATCH --account=stf-ckpt
#SBATCH --partition=ckpt

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

cat 32_fit_sharedq_jobs.sh  | parallel -j7
