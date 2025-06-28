#!/bin/bash

#SBATCH -o test_1.o%j

#SBATCH --time=72:00:00

#SBATCH -N 1 -n 15 --mem=150gb

#SBATCH --mail-user=cxm590@case.edu

#SBATCH --mail-type=all

#Copy to the scratch directory and change directory


cd $PWD


#Load MATLAB module

module load matlab


#RUN MATLAB script
  
matlab  -nodisplay -nosplash -softwareopengl -r  MyPredictivePriorCheck_1


echo Master process running on $HOSTNAME

echo Directory is `pwd`

NPROCS=$(( $SLURM_NPROCS*$SLURM_NNODES))

echo This job is allocated $NPROCS nodes

echo `date`

#Copy Back all the output files to the present working directory