#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=icelake
#PBS -A SlidingBootstrap
#PBS -N log2.0genCiBstMeanRFix


cd $PBS_O_WORKDIR #Verzeichnis von dem ich den Job aus starte

module load R/4.3.3-gcc



Rscript 2.0genCI_Cheap.R

