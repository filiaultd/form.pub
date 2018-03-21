#!/bin/bash

#PBS -P aquilegia
#PBS -N popool2
#PBS -j oe
#PBS -o output.popool2.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=62gb

### running popoolation2 for pools mapped to v3.1
### for individual formosa and pubescens pools

DATE=16Mar18
OUTDIR=/lustre/scratch/projects/aquilegia/016.popool2.all.pops.16March18/003.output
INDIR=/lustre/scratch/projects/aquilegia/016.popool2.all.pops.16March18/002.input
POOLDIR=/lustre/scratch/projects/aquilegia/001.common.reference.files/005.popoolation2
REF=/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta
OUT_PRE=no169.pools
SYNC_FILE=no169.pools.sync
MIN_C=4
MIN_COV=4
MAX_COV=200

module load Perl/5.20.0-goolf-1.4.10

cd $INDIR

# from previously-generated sync file, get AFDs

#perl $POOLDIR/snp-frequency-diff.pl --input $SYNC_FILE --output-prefix $OUT_PRE --min-count $MIN_C --min-coverage $MIN_COV --max-coverage $MAX_COV

# with AFDs, can do sliding window Fst


### large and small windows

#perl $POOLDIR/fst-sliding.pl --input  $OUT_PRE.sync --output $OUTDIR/$OUT_PRE.100kbwindow.50kbslide.fst --min-count $MIN_C --min-coverage $MIN_COV --max-coverage $MAX_COV --min-covered-fraction 0.0 --window-size 100000 --step-size 50000 --pool-size 40:47:22:26:48:49:24:68:37:35

#perl $POOLDIR/fst-sliding.pl --input  $OUT_PRE.sync --output $OUTDIR/$OUT_PRE.10kbwindow.5kbslide.fst --min-count $MIN_C --min-coverage $MIN_COV --max-coverage $MAX_COV --min-covered-fraction 0.0 --window-size 10000 --step-size 5000 --pool-size 40:47:22:26:48:49:24:68:37:35
 

### at every SNP:

perl $POOLDIR/fst-sliding.pl --input $OUT_PRE.sync --output $OUTDIR/$OUT_PRE.each.position.fst --suppress-noninformative --min-count $MIN_C --min-coverage $MIN_COV --max-coverage $MAX_COV --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 40:47:22:26:48:49:24:68:37:35



