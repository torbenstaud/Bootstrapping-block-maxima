#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=5GB:arch=skylake
#PBS -A SlidingBootstrap
#PBS -N z_genVarsNFix
#PBS -r y 
#PBS -J 1-168
#168 im echten

cd $PBS_O_WORKDIR #Verzeichnis von dem ich den Job aus starte

module load R/4.3.3-gcc

# Komm export R_LIBS=$R_LIBS:$HOME/R/libs
# env =^ Das sorgt dafür, dass alle 
#   Umgebungsvariablen des Moduls geprintet werden


Rscript 1generateTrueVarsMleBstRFix.R $PBS_ARRAY_INDEX

