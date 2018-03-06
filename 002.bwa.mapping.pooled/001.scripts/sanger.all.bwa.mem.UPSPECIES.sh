#!/bin/bash

#PBS -P aquilegia
#PBS -N bwa.NAMEID.map
#PBS -j oe
#PBS -o output.bwamap.JOBID.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=30gb

module load BWA/0.7.4-goolf-1.4.10
#module load SAMtools/0.1.18-goolf-1.4.10
#update samtools starting 24March15
module load SAMtools/0.1.19-goolf-1.4.10
module load GATK/2.5-2-Java-1.7.0_10
module load Java



######### IMPORTANT ####################################
# input is sanger encoded fastq files, unzipped.
# check the encoding of your files by running fastqc before proceeding
##################################################################

#cd /lustre/scratch/projects/aquilegia/000.gzipped.fastq

#REF=/lustre/scratch/projects/aquilegia/001.common.reference.files/009.aquilegia.assembly.20120630/002.new.bwa.index/Aquilegia_coerulea.main_genome.scaffolds.fasta
# new reference from 24 March 15
REF=/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta

upfile=UPSPECIES
#unzip_file=${upfile//.gz/}
unzip_file=$upfile
PREFIX=`echo $upfile | cut -d"." -f1`

DATE=24Mar15
GATKPATH=/lustre/scratch/projects/aquilegia/001.common.reference.files/002.GATK.source/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar

###################################
#### STEP 1: file prep and aligning
###################################

# 1.1 unzip file and move to correct directory
#gunzip $upfile
#echo $upfile
#echo $unzip_file
#echo $PREFIX
#mv $unzip_file /lustre/scratch/projects/aquilegia/005.mapping20May14
cd /lustre/scratch/projects/aquilegia/003.sanger.fastq.pools
##1.2 Start alignment

######################
# species stuff from JGI are interleaved fastq formatted and Illumina 1.5 encoded
######################
bwa mem -t 8 -p -M $REF $unzip_file > $PREFIX.$DATE.sam
samtools view -bhS -o $PREFIX.$DATE.bam $PREFIX.$DATE.sam
samtools sort $PREFIX.$DATE.bam $PREFIX.$DATE.sorted
samtools index $PREFIX.$DATE.sorted.bam
samtools rmdup $PREFIX.$DATE.sorted.bam $PREFIX.$DATE.dedup.bam
samtools index $PREFIX.$DATE.dedup.bam 

###############
#### NB: make sure the reference has a faidx and a dict file before using GATK!
#### Instructions:http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
################
 
#1.3 need to reheader bam
java -Xmx30g -jar /lustre/scratch/projects/aquilegia/001.common.reference.files/003.picardtools_1.53.source/AddOrReplaceReadGroups.jar INPUT=$PREFIX.$DATE.dedup.bam OUTPUT=$PREFIX.$DATE.dedup.reheadered.bam RGLB=$PREFIX RGPL=Illumina RGPU=$PREFIX RGSM=$PREFIX

#1.4 (optional) Quick coverage stuff for meeting tomorrow
samtools index $PREFIX.$DATE.dedup.reheadered.bam

java -Xmx30g -jar $GATKPATH \
   -R $REF \
   -T DepthOfCoverage \
   -o $PREFIX.$DATE.coverage.stats.txt \
   -I $PREFIX.$DATE.both.dedup.reheadered.bam

#1.5 remove intermediate files
if [ -e "$PREFIX.$DATE.dedup.reheadered.bam" ]
then
	rm $PREFIX.$DATE.sam
	rm $PREFIX.$DATE.sorted.bam
	rm $PREFIX.$DATE.sorted.bam.bai
	rm $PREFIX.$DATE.dedup.bam 
	rm $PREFIX.$DATE.dedup.bam.bai
	rm $unzip_file
fi
#files remaining are the original fastq, the original bam, and dedup.reheadered.bam (with .bai index)

################################################
## STEP TWO: clean up alignment to go into GATK
## GATK best practices
################################################

#2.1 clean bam file
java -Xmx30g -Dsnappy.disable=true -jar /lustre/scratch/projects/aquilegia/001.common.reference.files/003.picardtools_1.53.source/CleanSam.jar INPUT=$PREFIX.$DATE.dedup.reheadered.bam OUTPUT=$PREFIX.$DATE.clean.bam TMP_DIR=./999.tmp

#2.2 remove unmapped
samtools view -h -F 4 -b $PREFIX.$DATE.clean.bam > $PREFIX.$DATE.rm.unmap.bam

#2.3 resort
samtools sort $PREFIX.$DATE.rm.unmap.bam $PREFIX.$DATE.rm.unmap.sorted

#2.4 remove duplicates with picard
java -Xmx30g -Dsnappy.disable=true -jar /lustre/scratch/projects/aquilegia/001.common.reference.files/003.picardtools_1.53.source/MarkDuplicates.jar INPUT=$PREFIX.$DATE.rm.unmap.sorted.bam METRICS_FILE=$PREFIX.$DATE.dedup.metrics OUTPUT=$PREFIX.$DATE.pic.dedup.bam TMP_DIR=./999.tmp ASSUME_SORTED=TRUE

#2.5 fix mate pair issues
java -Xmx30g -Dsnappy.disable=true -jar /lustre/scratch/projects/aquilegia/001.common.reference.files/003.picardtools_1.53.source/FixMateInformation.jar INPUT=$PREFIX.$DATE.pic.dedup.bam OUTPUT=$PREFIX.$DATE.fixed.bam TMP_DIR=./999.tmp VALIDATION_STRINGENCY=LENIENT

#2.6 resort
samtools sort $PREFIX.$DATE.fixed.bam $PREFIX.$DATE.fixed.sorted

#2.7 index
samtools index $PREFIX.$DATE.fixed.sorted.bam

#2.8 removing unneeded intermediate files
if [ -e "$PREFIX.$DATE.fixed.sorted.bam.bai" ]
then
	rm $PREFIX.$DATE.dedup.reheadered.bam
        rm $PREFIX.$DATE.dedup.reheadered.bai
        rm $PREFIX.$DATE.clean.bam
	rm $PREFIX.$DATE.rm.unmap.bam
	rm $PREFIX.$DATE.rm.unmap.sorted.bam
	rm $PREFIX.$DATE.pic.dedup.bam
	rm $PREFIX.$DATE.fixed.bam
fi
# what remains are fixed.sorted.bam and fixed.sorted.bam.bai

###############################################
## STEP THREE: local realignments with GATK
###############################################

#3.1 id problem intervals
java -Xmx30g -jar $GATKPATH\
	-I $PREFIX.$DATE.fixed.sorted.bam\
	-R $REF\
	-T RealignerTargetCreator\
	-o $PREFIX.$DATE.intervals \

#3.2 actual realignments
java -Xmx30g -jar $GATKPATH\
	-I $PREFIX.$DATE.fixed.sorted.bam \
	-R $REF \
	-T IndelRealigner \
	-targetIntervals $PREFIX.$DATE.intervals \
	-o $PREFIX.$DATE.realigned.bam \

#3.3 index
samtools index $PREFIX.$DATE.realigned.bam 

#3.4 remove intermeidate files
if [ -e "$PREFIX.$DATE.realigned.bam.bai" ]
then
        rm $PREFIX.$DATE.fixed.sorted.bam
        rm $PREFIX.$DATE.fixed.sorted.bam.bai
        rm $PREFIX.$DATE.realigned.bai
fi


#########################################################################################################
## STEP FOUR: GATK SNP calling
## done in another script so can be done efficiently chromosome by chromosome.  See scripts in GATK folder in scripts.
#########################################################################################################





