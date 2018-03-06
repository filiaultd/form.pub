#!/bin/bash

#PBS -P aquilegia
#PBS -N to.pools
#PBS -j oe
#PBS -o output.bam.to.pools.log
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=6gb


#########################################################
#### taking individual mapped bam files and 
#### merging them by pool
#### a tedious proposition!
#### DLF 29Aug14
#########################################################

#module load SAMtools/0.1.18-goolf-1.4.10
#update samtools 26March15

module load SAMtools/0.1.19-goolf-1.4.10

cd /lustre/scratch/projects/aquilegia/003.sanger.fastq.pools

OUTPATH=/lustre/scratch/projects/aquilegia/003.sanger.fastq.pools/001.pool.bams


#pool1	Lundy Canyon Lower	form
#samtools merge $OUTPATH/pool1.bam 8281_sequence.24Mar15.realigned.bam 302GBAAXX_1_8281_20091109_seq.24Mar15.realigned.bam

#pool2	Rock Creek		form
#samtools merge $OUTPATH/pool2.bam 8282_sequence.24Mar15.realigned.bam 302GBAAXX_2_8282_20091109_seq.24Mar15.realigned.bam

#pool3	Gable Lakes		form
#samtools merge $OUTPATH/pool3.bam 8434_sequence.24Mar15.realigned.bam s_1_sequence.24Mar15.realigned.bam

#pool4	North Lake Stream	form
#samtools merge $OUTPATH/pool4.bam 8283_sequence.24Mar15.realigned.bam 302GBAAXX_3_8283_20091109_seq.24Mar15.realigned.bam

#pool5	Bishop Pass Trail	form
#samtools merge $OUTPATH/pool5.bam 8284_sequence.24Mar15.realigned.bam 302GBAAXX_4_8284_20091109_seq.24Mar15.realigned.bam

#pool6	Golden Trout 1		form
#samtools merge $OUTPATH/pool6.bam 8285_sequence.24Mar15.realigned.bam 302GBAAXX_5_8285_20091109_seq.24Mar15.realigned.bam

#pool7	Big Horn Lake		pub
#samtools merge $OUTPATH/pool7.bam 8435_sequence.24Mar15.realigned.bam s_2_sequence.24Mar15.realigned.bam

#pool9	Morgan Pass		pub
samtools merge $OUTPATH/pool9.bam 8288_sequence.24Mar15.realigned.bam 302GBAAXX_6_8288_20091109_seq.24Mar15.realigned.bam s_6_sequence.24Mar15.realigned.bam

#pool10	Piute Pass		pub
#samtools merge $OUTPATH/pool10.bam 8289_sequence.24Mar15.realigned.bam 302GBAAXX_7_8289_20091109_seq.24Mar15.realigned.bam

#pool11	Bishop Pass		pub
#samtools merge $OUTPATH/pool11.bam 8290_sequence.24Mar15.realigned.bam 302GBAAXX_8_8290_20091109_seq.24Mar15.realigned.bam

#pool12	Golden Trout 2		form
#cp 81KHFABXX_7_31-08-2011.24Mar15.realigned.bam $OUTPATH/pool12.bam

#pool13	Lamark Col		pub
#cp 81KHFABXX_6_31-08-2011.24Mar15.realigned.bam $OUTPATH/pool13.bam

#pool14	Mono Pass		pub
cp C0HWGACXX_1_20120713B_20120724.24Mar15.realigned.bam $OUTPATH/pool14.bam
