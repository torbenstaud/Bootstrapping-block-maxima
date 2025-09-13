#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=12GB:arch=icelake
#PBS -A SlidingBootstrap
#PBS -N generateMean
#PBS -r y 
#PBS -J 1-20


cd $PBS_O_WORKDIR 

module load R/4.3.3-gcc


Rscript 1genDataHpcMean.R $PBS_ARRAY_INDEX

