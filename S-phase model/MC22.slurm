#!/bin/bash

#SBATCH -o test.o%j

#SBATCH --time=100:00:00

#SBATCH -N 1 -n 10 --mem=50gb

#SBATCH --mail-user=cxm590@case.edu

#SBATCH --mail-type=all

#Copy to the scratch directory and change directory


cd $PWD


#Load MATLAB module

module load matlab


#RUN MATLAB script
  
matlab  -nodisplay -nosplash -softwareopengl -r S_phase_MCMC_main_treated


echo Master process running on $HOSTNAME

echo Directory is `pwd`

NPROCS=$(( $SLURM_NPROCS*$SLURM_NNODES))

echo This job is allocated $NPROCS nodes

echo `date`

#Copy Back all the output files to the present working directory
