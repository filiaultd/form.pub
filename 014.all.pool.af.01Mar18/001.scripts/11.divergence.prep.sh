#!/bin/bash

#PBS -P aquilegia
#PBS -N sync.parse
#PBS -j oe
#PBS -o output.sync.parse.JOBID.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=60gb


scriptPath=/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts
syncPath=/lustre/scratch/projects/aquilegia/001.pool.bams
syncFile=no169.pools.sync
#syncFile=test.sync
outPath=/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input
outFile=no169.pools.poly.sync
minCov=4
minPass=4

cd $scriptPath

python 10.divergence.prep.py $syncPath/$syncFile $outPath/$outFile $minCov $minPass
 
