#!/bin/bash

#PBS -P aquilegia
#PBS -N polar.prep
#PBS -j oe
#PBS -o output.polar.prep.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=32GB:local_disk=30gb

# shell script for filtering vcf and getting genotype table
# for Semi data (mapping is from the AQ genomes paper)
# overall goal is to polarize pooled data vs. Semi
# and look for sharing patterns of SNPs with high AFD
# DLF 07May18

module load Java/1.7.0_10
# set some variables here
TEMPPATH=$LOCAL_DISK
REF=/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta
DATE=07May18
gVCF=SRR4.20July15.combined.g.vcf
VCFUP=SRR4.$DATE.vcf
VCFOUT=SRR4.$DATE.pool.subset.vcf
POOLSNPS_FILE=no169.pools.poly.sync.biallelic.out
POOLSNPS=pool.biallelic.SNPS.table
GATKPATH=/lustre/scratch/projects/aquilegia/001.common.reference.files/002.GATK.source/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar
OUTPATH=/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input/
VCFPATH=$OUTPATH

cd $VCFPATH

###########################
### get ready
### re-genotype Semi (SRR4) from g.vcf
############################

#java -Xmx30G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
#        -R $REF\
#        -T GenotypeGVCFs\
#        -V $gVCF\
#        --includeNonVariantSites \
#        -o $VCFUP

####################################
### first step, mask positions not in 
### pool biallelic SNPs
######################################

#first make table of sites to extract from Semi snps
#format is first row "HEADER"
#next rows, one site per row with format chr:pos
awk '{print $1 ":" $2}' $POOLSNPS_FILE > tmp.txt
sed -i '1s/^/HEADER\n/' tmp.txt
mv tmp.txt $POOLSNPS

java -Xmx30G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
       -R $REF\
       -T VariantFiltration\
	-V $VCFUP\
	-o $VCFOUT\
	-mask:TABLE $POOLSNPS\
	-filterNotInMask \
	-maskName NO.POOL

#####################################
### final step
### output to table
### shouldn't include filtered sites from filter above
######################################


java -Xmx30G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
        -R $REF\
        -T VariantsToTable\
        -V $VCFOUT\
        -F CHROM -F POS -F REF -F ALT -F FILTER -F NO-CALL\
	-GF GT \
        --allowMissingData \
        -o SRR.pool.snps.genotype.table.txt
