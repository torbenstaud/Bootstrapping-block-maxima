#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=32GB:arch=skylake
#PBS -A SlidingBootstrap
#PBS -N bstMleEstData
#PBS -r y 
#PBS -J 1-24

cd $PBS_O_WORKDIR 

module load R/4.3.3-gcc

Rscript 1genDataHpcRFix.R $PBS_ARRAY_INDEX

