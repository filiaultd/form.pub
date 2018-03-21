#!/bin/bash

#PBS -P aquilegia
#PBS -N popool2gene
#PBS -j oe
#PBS -o output.popool2gene.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=30gb

## doing popoolation2 for mapping with v3 genome using new v3.1 annotation

module load Java/1.7.0_21
module load Perl/5.20.0-goolf-1.4.10

### set variables here

DATE=15MAR18
POOL2_DIR=/lustre/scratch/projects/aquilegia/001.common.reference.files/005.popoolation2
INDIR=/lustre/scratch/projects/aquilegia/016.popool2.all.pops.16March18/002.input
SYNC_FILE=no169.pools.sync
GTF_FILE=/lustre/scratch/projects/aquilegia/001.common.reference.files/018.Acoerulea_v3.1_annotation_20150930/Acoerulea_396_v3.1_genes_only.gtf
OUTPUT_DIR=/lustre/scratch/projects/aquilegia/016.popool2.all.pops.16March18/003.output
OUT_PREFIX=no169.genes.$DATE
SYNC_GENE_FILE=$OUT_PREFIX.sync
FST_FILE=$OUT_PREFIX.fst
MIN_C=4
MIN_COV=4
MAX_COV=200



cd $INDIR


## STEP ONE: make sync file.  (Alredy have mpileup file for this mapping, AND sync file, so can skip this step for now, which is good because it is slow)

## STEP TWO: make gene-wise sync file
#perl $POOL2_DIR/create-genewise-sync.pl --input $SYNC_FILE --gtf $GTF_FILE --output $SYNC_GENE_FILE

## STEP THREE: get frequency diffs
#perl $POOL2_DIR/snp-frequency-diff.pl --input $SYNC_GENE_FILE --output-prefix $OUT_PREFIX --min-count $MIN_C --min-coverage $MIN_COV --max-coverage $MAX_COV

## STEP FOUR: calculate FST
perl $POOL2_DIR/fst-sliding.pl --input $SYNC_GENE_FILE  --output $OUTPUT_DIR/$FST_FILE --min-count $MIN_C --min-coverage $MIN_COV --max-coverage $MAX_COV --min-covered-fraction 0.0 --window-size 100000 --step-size 100000 --pool-size 40:47:22:26:48:49:24:68:37:35

