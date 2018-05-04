#!/bin/bash

#PBS -P aquilegia
#PBS -N parse.afs
#PBS -j oe
#PBS -o output.get.afs.JOBID.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=100gb



scriptPath=/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts

cd $scriptPath

module load R

Rscript 14.parse.afd.R
 
