#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=32GB:arch=skylake
#PBS -A SlidingBootstrap
#PBS -N bstMleBstData
#PBS -r y 
#PBS -J 1-192
#J war vorher 1-2; 

cd $PBS_O_WORKDIR 

module load R/4.3.3-gcc



Rscript 1genDataHpcMleBstNFix.R $PBS_ARRAY_INDEX

