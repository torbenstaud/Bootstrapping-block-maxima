#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=10GB:arch=icelake
#PBS -A SlidingBootstrap
#PBS -N bstMLE


cd $PBS_O_WORKDIR #Verzeichnis von dem ich den Job aus starte

module load R/4.3.3-gcc

# Komm export R_LIBS=$R_LIBS:$HOME/R/libs
# env =^ Das sorgt dafür, dass alle 
#   Umgebungsvariablen des Moduls geprintet werden


Rscript 2.1procDataHpcMean.R

