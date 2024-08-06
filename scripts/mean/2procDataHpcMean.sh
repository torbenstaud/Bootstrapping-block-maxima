#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=icelake
#PBS -A SlidingBootstrap
#PBS -N procMean
#PBS -r y 
#PBS -J 1-20

cd $PBS_O_WORKDIR #Verzeichnis von dem ich den Job aus starte

module load R/4.3.3-gcc

# Komm export R_LIBS=$R_LIBS:$HOME/R/libs
# env =^ Das sorgt dafür, dass alle 
#   Umgebungsvariablen des Moduls geprintet werden


Rscript 2procDataHpcMean.R $PBS_ARRAY_INDEX

