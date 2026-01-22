#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=u1-qlm
#SBATCH --output=slurm_logs/slurm-%j.out

WORKDIR=/space2/kanwar/u1-qlm-heavy-hex
srun "$@"
