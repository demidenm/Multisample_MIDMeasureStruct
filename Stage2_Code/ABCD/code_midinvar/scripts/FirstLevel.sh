#!/bin/bash -l
#SBATCH -J first_projinv
#SBATCH --array=0-779 # jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6G
#SBATCH -t 00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}.edu
#SBATCH -p msismall,amdsmall
#SBATCH -o batch_logs/%x_%A_%a.out
#SBATCH -e batch_logs/%x_%A_%a.err
#SBATCH -A ${PROFILE}

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env
module load fsl

ID=${SLURM_ARRAY_TASK_ID}

bash ./batch_run/first${ID}
