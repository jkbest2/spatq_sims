#! /bin/bash

qdevscalingsim_jid=$(sbatch 21_runsims_qdevscaling.slurm)
qdevscalingfit_jid=$(sbatch --dependency=afterok:$qdevscalingsim_jid 31_fit_qdevscaling_parallel_ckpt.slurm)

sharedqsim_jid=$(sbatch 21_runsims_sharedq.slurm)
sharedqfit_jid=$(sbatch --dependency=afterok:$sharedqsim_jid 31_fit_sharedq_parallel_ckpt.slurm)

prefintensitysim_jid=$(sbatch 21_runsims_prefintensity.slurm)
prefintensityfit_jid=$(sbatch --dependency=afterok:$prefintensitysim_jid 31_fit_prefintensity_parallel_ckpt.slurm)

sbatch --dependency=afterok:$qdevscalingfit_jid:$sharedqfit_jid:$prefintensityfit_jid 40_readfits.slurm

# show dependencies in squeue output:
squeue -u $USER -o "%.8A %.4C %.10m %.20E"
