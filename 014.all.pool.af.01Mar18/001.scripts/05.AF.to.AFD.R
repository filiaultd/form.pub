### from allele counts, get allele frequency difference and p-values
### DLF 15MAR18

InDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output"
FunFile="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/04.AF.helper.fxns.R"
AFDatFile="biallele.freq.Rdata"
BamDir="/lustre/scratch/projects/aquilegia/001.pool.bams"
BamList="no169.pool.nodups.bam.list"
OutFile="biallele.afs.Rdata"


setwd(InDataDir)

source(FunFile)

load(AFDatFile)


#############################################
### step three: AFD between species pools ###
#############################################
pools <- read.table(paste(BamDir, BamList, sep="/"),stringsAsFactors=FALSE)
pools <- matrix(pools[,1],nrow=10, ncol=1)
pools <- gsub(".bam.no.dups","", pools)

biallele.freq <- as.data.frame(biallele.freq, stringsAsFactors=FALSE)
colnames(biallele.freq) <- c("chr","pos","ref","alleles",pools)
pub.pools <- paste("pool",c(7,10,11,13,14),sep="")
form.pools <- paste("pool",c(2,3,4,5,12), sep="")
sp.index <- cbind(c(pub.pools,form.pools), c(rep("p",length(pub.pools)), rep("f", length(form.pools))))
sp.index <- as.data.frame(sp.index, stringsAsFactors=FALSE)
colnames(sp.index) <- c("pool","species")

#test <- apply(biallele.freq[1:100,],1,get.af)

biallele.afs <- apply(biallele.freq, 1, get.af)
save(biallele.afs, file=OutFile)


##############################
### write chr and pos of biallele.freq
### to use as a filter for Semi SNPs vcf in GATK
#########################################

options(scipen = 999)
bi.pos <- paste(biallele.freq[,1], as.numeric(biallele.freq[,2]), sep=":")
bi.pos <- c("HEADER", bi.pos)
write.table(bi.pos, file="pool.biallelic.SNPS.table", quote=FALSE,col.names=FALSE, row.names=FALSE)



################## function development below this line ############################


get.af <- function(af.line){
	alleles <- unlist(af.line[5:14])
	alleles <- strsplit(alleles, "/")
	alleles <- do.call(rbind, alleles)
	class(alleles) <- "numeric"
	colnames(alleles) <- unlist(strsplit(af.line[4],"/"))
	alleles <- as.data.frame(alleles)
	alleles$species <- sp.index[match(rownames(alleles),sp.index$pool),2]
	
	tt1 <- try(t.test(alleles[,1]~alleles[,3], data=alleles))
	if(length(tt1)==1){
		af1 <- aggregate(alleles[,1]~alleles[,3], data=alleles, mean)
		t1.sum <- c(af1[,2],NA)
	}else{	
		t1 <- t.test(alleles[,1]~alleles[,3], data=alleles)
		t1.sum <- c(t1$estimate, t1$p.value)
		}
	names(t1.sum) <- c("f1","p1","pval1")
	
	tt2 <- try(t.test(alleles[,2]~alleles[,3], data=alleles))
	if(length(tt2)==1){
		af2 <- aggregate(alleles[,2]~alleles[,3], data=alleles, mean)
		t2.sum <- c(af2[,2], NA)
	}else{
		t2 <- t.test(alleles[,2]~alleles[,3], data=alleles)
		t2.sum <- c(t2$estimate, t2$p.value)
		}
	names(t2.sum) <- c("f2","p2","pval2")
	
	out.dat <- unlist(c(t1.sum, t2.sum))
	return(out.dat)
}
