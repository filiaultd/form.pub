#!/bin/bash

#PBS -P aquilegia
#PBS -N bam.to.mpileup
#PBS -j oe
#PBS -o /lustre/scratch/projects/aquilegia/001.pool.bams
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=16:mem=60GB:local_disk=60gb



### shell script for combining each pool bam file into an mpileup and then sync file
### prepwork for popoolation 2

cd /lustre/scratch/projects/aquilegia/001.pool.bams

module load SAMtools/1.6-foss-2017a
module load Java

BAM_LIST=no169.pool.bams.list

### pools 1,6,9 have very low coverage - much lower than others!
### keep 4 pools of each species
### see previous analyses of why no 1,6,9

#### first step is to remove ambiguous reads

while read UPBAM
do
    echo $UPBAM
    samtools view -q20 -bS $UPBAM | samtools sort - -o $UPBAM.no.dups
done < $BAM_LIST

#### make these into mpileup

REF_FASTA=/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta
OUT_mpileup=no169.pools.mpileup
ls *no.dups > no169.pool.nodups.bam.list
ND_BAM_LIST=no169.pool.nodups.bam.list

samtools mpileup -f $REF_FASTA -b $ND_BAM_LIST -o $OUT_mpileup


##### generate sync file for popoolation2

PP2_DIR=/lustre/scratch/projects/aquilegia/001.common.reference.files/005.popoolation2
OUT_sync=no169.pools.sync

#perl $PP2_DIR/mileup2sync.pl --fastq-type sanger --min-qual 20 --input $OUT_mpileup --output $OUT_sync
java -ea -Xmx8g -jar $PP2_DIR/mpileup2sync.jar --input $OUT_mpileup --output $OUT_sync --fastq-type sanger --min-qual 20 --threads 8
