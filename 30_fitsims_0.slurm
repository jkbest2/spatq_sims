#!/bin/bash

## Job name
#SBATCH --job-name=fit_sims

## Allocations definition
#SBATCH --account=stf
#SBATCH --partition=stf

## Nodes
#SBATCH --nodes=1

## Tasks
#SBATCH --tasks=1

## CPUs per task
#SBATCH --cpus-per-task=8

## Walltime
#SBATCH --time=16:00:00

## Memory per node
#SBATCH --mem=40G

## Set the working directory for this job.
## Directory must exist before submitting job.
#SBATCH --chdir=/gscratch/stf/jkbest/spatq_sims

## Use email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jkbest@uw.edu

Rscript 30_fitsims.R 20 25
