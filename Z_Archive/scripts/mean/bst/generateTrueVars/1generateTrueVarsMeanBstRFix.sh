#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=skylake
#PBS -A SlidingBootstrap
#PBS -N z_genVarsNFix
#PBS -r y 
#PBS -J 1-280


cd $PBS_O_WORKDIR

module load R/4.3.3-gcc



Rscript 1generateTrueVarsMeanBstRFix.R $PBS_ARRAY_INDEX

