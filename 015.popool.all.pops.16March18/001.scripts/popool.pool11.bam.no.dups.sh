#!/bin/bash

#PBS -P aquilegia
#PBS -N pop.pool11.bam.no.dups
#PBS -j oe
#PBS -o output.popool.pool11.bam.no.dups
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=4:mem=30gb

module load SAMtools/0.1.19-goolf-1.4.10

### doing popoolation calculations of pi and theta
### can also do tajima's D, but not yet run
### for each individual pool
### two windows, individual SNP, by gene
### DLF 16March18

UPDIR=/lustre/scratch/projects/aquilegia/015.popool.all.pops.16March18
DATE=16MAR18

BAMDIR=/lustre/scratch/projects/aquilegia/001.pool.bams
BAMUP=pool11.bam.no.dups
MPILEUPFILE=$BAMUP.mpileup
PI_OUT_FILE=$BAMUP.pi.out
THETA_OUT_FILE=$BAMUP.theta.out
DOWN_MPILEUP=$BAMUP.downsampled.mpileup
TAJD_OUT_FILE=$BAMUP.tajimasd.out
POOL_SIZE=47
POOL_MEDIAN=9
INDIR=$UPDIR/002.input
OUTDIR=$UPDIR/003.output
POPOOLDIR=/lustre/scratch/projects/aquilegia/001.common.reference.files/016.popoolation
REF=/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta
GTFFILE=/lustre/scratch/projects/aquilegia/001.common.reference.files/018.Acoerulea_v3.1_annotation_20150930/Acoerulea_396_v3.1_genes_only.gtf


cd $UPDIR

# 1. generate mpileup
samtools mpileup -f $REF $BAMDIR/$BAMUP > $INDIR/$MPILEUPFILE

cd $INDIR

# 2. run popoolation pi in sliding window
# three window sizes

## 100kb with 50kb slide
perl $POPOOLDIR/Variance-sliding.pl --input $MPILEUPFILE --output $OUTDIR/$PI_OUT_FILE.100kb_50kb --measure pi --window-size 100000 --step-size 50000 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --pool-size $POOL_SIZE --fastq-type sanger

## 10kb with 5kb slide
perl $POPOOLDIR/Variance-sliding.pl --input $MPILEUPFILE --output $OUTDIR/$PI_OUT_FILE.10kb_5kb --measure pi --window-size 10000 --step-size 5000 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --pool-size $POOL_SIZE --fastq-type sanger

## individual SNPs
perl $POPOOLDIR/Variance-sliding.pl --input $MPILEUPFILE --output $OUTDIR/$PI_OUT_FILE.1b_1b --measure pi --window-size 1 --step-size 1 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --pool-size $POOL_SIZE --fastq-type sanger



# 3. run theta in sliding window
# three window sizes

## 100kb with 50kb slide
perl $POPOOLDIR/Variance-sliding.pl --input $MPILEUPFILE --output $OUTDIR/$PI_OUT_FILE.100kb_50kb --measure theta --window-size 100000 --step-size 50000 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --pool-size $POOL_SIZE --fastq-type sanger

## 10kb with 5kb slide
perl $POPOOLDIR/Variance-sliding.pl --input $MPILEUPFILE --output $OUTDIR/$PI_OUT_FILE.10kb_5kb --measure theta --window-size 10000 --step-size 5000 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --pool-size $POOL_SIZE --fastq-type sanger

## individual SNPs
perl $POPOOLDIR/Variance-sliding.pl --input $MPILEUPFILE --output $OUTDIR/$PI_OUT_FILE.1b_1b --measure theta --window-size 1 --step-size 1 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --pool-size $POOL_SIZE --fastq-type sanger


#4. run pi per gene
perl $POPOOLDIR/Variance-at-position.pl --measure pi --pileup $MPILEUPFILE --gtf $GTFFILE --output $OUTDIR/gene.$PI_OUT_FILE --pool-size $POOL_SIZE --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --fastq-type sanger

# 5. run theta per gene
perl $POPOOLDIR/Variance-at-position.pl --measure theta --pileup $MPILEUPFILE --gtf $GTFFILE --output $OUTDIR/gene.$THETA_OUT_FILE --pool-size $POOL_SIZE --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --fastq-type sanger

########## may want to think about the max coverage thing a bit
######### based on coverage per pool?
######### 70 is probably very high

# 4. run Tajima's D in sliding window

# 4a. first step is to downsample so even coverage - use median as calculated from BAM file

#perl $POPOOLDIR/basic-pipeline/subsample-pileup.pl --input $MPILEUPFILE --output $DOWN_MPILEUP --target-coverage $POOL_MEDIAN --max-coverage 70 --fastq-type sanger --method withoutreplace --min-qual 20 

# 4b. second step is to calculate by gene 

#perl $POPOOLDIR/Variance-at-position.pl --measure D --pileup $DOWN_MPILEUP --gtf $GTFFILE --output $OUTDIR/gene.$TAJD_OUT_FILE --pool-size $POOL_SIZE --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20 --fastq-type sanger
