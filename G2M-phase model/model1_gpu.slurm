#!/bin/bash

#SBATCH -o test.o%j
#SBATCH --time=300:00:00
#SBATCH -N 1 -n 1 --mem=150gb
#SBATCH -p gpu  --gres=gpu:1                  # Request 1 GPU
#SBATCH --mail-user=cxm590@case.edu
#SBATCH --mail-type=ALL

echo Master process running on $HOSTNAME
echo Directory is `pwd`
echo `date`

# Load modules
module load matlab

# Run MATLAB
matlab -nodisplay -nosplash -softwareopengl -r "MCMC_main_treated"

# Copy results (optional)
# cp -r $SLURM_TMPDIR/* $SLURM_SUBMIT_DIR/
