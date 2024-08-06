#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=skylake
#PBS -A SlidingBootstrap
#PBS -N bstMleEstData
#PBS -r y 
#PBS -J 1-168

cd $PBS_O_WORKDIR 

module load R/4.3.3-gcc


Rscript 1genDataHpcMleBstRFix.R $PBS_ARRAY_INDEX

