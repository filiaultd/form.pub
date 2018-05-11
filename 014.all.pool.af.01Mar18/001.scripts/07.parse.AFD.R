### parsing allele frequency differences
### DLF 15Mar18


#### this seems to be obsolete??

InDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output"
FunFile="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/04.AF.helper.fxns.R"
BamDir="/lustre/scratch/projects/aquilegia/001.pool.bams"
BamList="no169.pool.nodups.bam.list"
InDataFile="biallele.afs.Rdata"

setwd(InDataDir)

source(FunFile)

load(InDataFile)  #biallele.afs

biallele.afs <- t(biallele.afs)

afd1 <- biallele.afs[,1]-biallele.afs[,2]
afd2 <- biallele.afs[,4]-biallele.afs[,5]
