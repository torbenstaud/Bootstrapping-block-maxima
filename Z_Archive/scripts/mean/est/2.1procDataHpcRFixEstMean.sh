#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=icelake
#PBS -A SlidingBootstrap
#PBS -N 2.1MeanEstProc



cd $PBS_O_WORKDIR

module load R/4.3.3-gcc


Rscript 2.1procDataHpcRFixEstMean.R

