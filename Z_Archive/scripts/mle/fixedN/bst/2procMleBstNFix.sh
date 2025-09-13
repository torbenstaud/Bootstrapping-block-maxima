#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=icelake
#PBS -A SlidingBootstrap
#PBS -N bstMle2RNFixProc


cd $PBS_O_WORKDIR 

module load R/4.3.3-gcc



Rscript 2procMleBstNFix.R

