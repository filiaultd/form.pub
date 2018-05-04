#!/bin/bash

#PBS -P aquilegia
#PBS -N AFtoAFD
#PBS -j oe
#PBS -o AF.to.AFD.JOBID.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=150gb


scriptPath=/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts
scriptFile=05.AF.to.AFD.R

cd $scriptPath

module load R

Rscript $scriptPath/$scriptFile

 
