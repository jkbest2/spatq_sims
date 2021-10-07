#! /bin/bash

qdevscalingsim_jid=$(sbatch --parsable 21_runsims_qdevscaling.slurm)
qdevscalingfit_jid=$(sbatch --parsable --dependency=afterok:$qdevscalingsim_jid 31_fit_qdevscaling_ckpt.slurm)
sbatch --dependency=singleton 31_fit_qdevscaling_ckpt.slurm
sbatch --dependency=singleton 31_fit_qdevscaling_ckpt.slurm
sbatch --dependency=singleton 31_fit_qdevscaling_ckpt.slurm
sbatch --dependency=singleton 31_fit_qdevscaling_ckpt.slurm

# sharedqsim_jid=$(sbatch --parsable 22_runsims_sharedq.slurm)
# sharedqfit_jid=$(sbatch --parsable --dependency=afterok:$sharedqsim_jid 32_fit_sharedq_ckpt.slurm)
# sbatch --dependency=singleton 32_fit_sharedq_ckpt.slurm
# sbatch --dependency=singleton 32_fit_sharedq_ckpt.slurm
# sbatch --dependency=singleton 32_fit_sharedq_ckpt.slurm
# sbatch --dependency=singleton 32_fit_sharedq_ckpt.slurm

prefintensitysim_jid=$(sbatch --parsable 23_runsims_prefintensity.slurm)
prefintensityfit_jid=$(sbatch --parsable --dependency=afterok:$prefintensitysim_jid 33_fit_prefintensity_ckpt.slurm)
sbatch --dependency=singleton 33_fit_prefintensity_ckpt.slurm
sbatch --dependency=singleton 33_fit_prefintensity_ckpt.slurm
sbatch --dependency=singleton 33_fit_prefintensity_ckpt.slurm
sbatch --dependency=singleton 33_fit_prefintensity_ckpt.slurm

densdepqsim_jid=$(sbatch --parsable 24_runsims_densdepq.slurm)
densdepqfit_jid=$(sbatch --parsable --dependency=afterok:$densdepqsim_jid 34_fit_densdepq_ckpt.slurm)
sbatch --dependency=singleton 34_fit_densdepq_ckpt.slurm
sbatch --dependency=singleton 34_fit_densdepq_ckpt.slurm
sbatch --dependency=singleton 34_fit_densdepq_ckpt.slurm
sbatch --dependency=singleton 34_fit_densdepq_ckpt.slurm

# sbatch --dependency=afterok:$qdevscalingfit_jid:$sharedqfit_jid:$prefintensityfit_jid:$densdepqfit_jid 40_readfits.slurm

# show dependencies in squeue output:
squeue -u $USER -o "%.8A %.4C %.10m %.20E"
