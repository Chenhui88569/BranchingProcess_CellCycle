#!/bin/bash

#SBATCH -o test.o%j

#SBATCH --time=100:00:00

#SBATCH  -n 30 --mem=200gb

#SBATCH --mail-user=cxm590@case.edu

#SBATCH --mail-type=all

#Copy to the scratch directory and change directory


cd $PWD


#Load MATLAB module

module load  matlab/R2021b


#RUN MATLAB script
  
matlab  -nodisplay -nosplash -softwareopengl -r  CredibleIntervalCompute


echo Master process running on $HOSTNAME

echo Directory is `pwd`

NPROCS=$(( $SLURM_NPROCS*$SLURM_NNODES))

echo This job is allocated $NPROCS nodes

echo `date`

#Copy Back all the output files to the present working directory
