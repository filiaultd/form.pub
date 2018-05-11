#### script to polarize SNPs and calculate allele frequencies of derived alleles
#### prep work to id bialleleic SNPs is done in script 10 and 11 in this directory
#### prep work to generate Semi SNPs is from script 24
#### DLF 08May18


SemiGtTable="SRR.pool.snps.genotype.table.txt"
InDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input"
InDataFile="no169.pools.poly.sync.biallelic.out"
OutDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output"
min.MAC=1
FunFile="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/04.AF.helper.fxns.R"
BamDir="/lustre/scratch/projects/aquilegia/001.pool.bams"
BamList="no169.pool.nodups.bam.list"


setwd(InDataDir)
source(FunFile)


#############################################
### get pop frequencies of biallelic SNPs
#################################################

bi.snps <- read.table(InDataFile, stringsAsFactors=FALSE)  ### slightly different from Semi - about 80K more positions, but I'm not too worried about this.
bi.snps <- bi.snps[grep("Chr", bi.snps[,1]),]

pools <- read.table(paste(BamDir, BamList, sep="/"),stringsAsFactors=FALSE)
pools <- matrix(pools[,1],nrow=10, ncol=1)
pools <- gsub(".bam.no.dups","", pools)
pub.pools <- paste("pool",c(7,10,11,13,14),sep="")
form.pools <- paste("pool",c(2,3,4,5,12), sep="")
sp.index <- cbind(c(pub.pools,form.pools), c(rep("p",length(pub.pools)), rep("f", length(form.pools))))
sp.index <- as.data.frame(sp.index, stringsAsFactors=FALSE)
colnames(sp.index) <- c("pool","species")
colnames(bi.snps) <- c("CHROM","POS","REF",pools)


#bi.dat.freq <- apply(bi.snps[1:1000,],1, function(x){parse.line(x,min.MAC)})
bi.dat.freq <- apply(bi.snps,1, function(x){parse.line(x,min.MAC)})
bi.dat.freq <- t(bi.dat.freq)
bi.dat.freq <- as.data.frame(bi.dat.freq, stringsAsFactors=FALSE)
class(bi.dat.freq$POS) <- "numeric"

save(bi.dat.freq, file=paste(OutDataDir,"bi.dat.freq.Rdata", sep="/"))


#####################
### load Semi snps
######################

semi.snps <- read.table(SemiGtTable, stringsAsFactors=FALSE, header=TRUE) #24435765 lines

### 5005716 uncalled in Semi, 19430049 called positions
#table(semi.snps$SRR4.GT) -> ugh
#ugh <- ugh[order(ugh, decreasing=TRUE)]
#blah <- sapply(names(ugh), nchar)
#bi.names <- ugh[blah==3]
#sum(bi.names)-(5005716+1)  ### number ./. or N/N
#### 19351573 called and bialleleic SNP (or nonpolymorphic)

### get positions where Semi is called and het
hom <- c("A/A","T/T","C/C","G/G")
g.semi.snps <- semi.snps[semi.snps$SRR4.GT%in%hom,]  
g.semi.snps <- g.semi.snps[grep("Chr", g.semi.snps[,1]),]  #18392681positions

########################
### get derived allele frequency for all SNPs
########################## 

all.dat <- merge(bi.dat.freq, g.semi.snps, all.x=TRUE)
all.dat.semi <- all.dat[is.na(all.dat$SRR4.GT)==FALSE,]
all.dat.semi$sa <- sapply(all.dat.semi$SRR4.GT, function(x) substr(x,1,1))


derived.af <- apply(all.dat.semi,1,get.derived.af)
derived.af <- t(derived.af)
colnames(derived.af) <- pools
save(derived.af, file="derived.af.Rdata")


p.pools <- sp.index[sp.index$species=="p",1]
p.daf <- derived.af[,colnames(derived.af)%in%p.pools]
p.mean <- apply(p.daf,1, function(x) mean(x, na.rm=TRUE))
f.pools <- sp.index[sp.index$species=="f",1]
f.daf <- derived.af[,colnames(derived.af)%in%f.pools]
f.mean <- apply(f.daf,1, function(x) mean(x, na.rm=TRUE))
both.mean.af <- cbind(p.mean, f.mean)
# 18150806 SNPs total polarized


##################################
### plot joint SFS of pub, form
#####################################

pdf("../003.output/derived.joint.SFS.pdf")
both.mean.af <- cbind(p.mean, f.mean)
smoothScatter(both.mean.af[,1], both.mean.af[,2], xlab="formosa mean allele frequency", ylab="pubescens mean allele frequency")
dev.off()

######################################
### to simplify, do this in bins
########################################

u <- 0.05 #bin size
bin.af <- apply(both.mean.af, 2, function(x){floor(x/u)*u}) # rounds down to nearest bin.size
colnames(bin.af) <- c("pub","form")
bin.af <- as.data.frame(bin.af)
bin.sfs <- with(bin.af, table(pub,form))

write.table(bin.sfs, file="../003.output/SFS.derived.form.pub.0.05bin.txt")
pdf(file="../003.output/SFS.derived.form.pub.0.05bin.pdf", width=6, height=6)
image(log(bin.sfs), xlab="pubescens derived allele freq", ylab="formosa derived allele freq", col=colorRampPalette(c("lightsteelblue1","steelblue4"))(10))
dev.off()

u <- 0.1 #bin size
bin.af <- apply(both.mean.af, 2, function(x){floor(x/u)*u}) # rounds down to nearest bin.size
colnames(bin.af) <- c("pub","form")
bin.af <- as.data.frame(bin.af)
bin.sfs <- with(bin.af, table(pub,form))

write.table(bin.sfs, file="../003.output/SFS.derived.form.pub.0.1bin.txt")
pdf(file="../003.output/SFS.derived.form.pub.0.1bin.pdf", width=6, height=6)
image(log(bin.sfs), xlab="pubescens derived allele freq", ylab="formosa derived allele freq",col=colorRampPalette(c("lightsteelblue1","steelblue4"))(10))
dev.off()

###############################
### get some quick ideas of AFD numbers
#######################################
fp.diff <- f.mean-p.mean

fixed <- fp.diff[abs(fp.diff)==1]
fixed <- fixed[is.na(fixed)==FALSE]
table(fixed)
#fixed
# -1   1 
# 98 251 
# so slightly more fixed in formosa than pub, too.

diff9 <- fp.diff[abs(fp.diff)>0.9]
diff9 <- diff9[is.na(diff9)==FALSE] #1463


load("25.examine.polarized.allele.freqs.Rdata")

#####################################
### do some plotting locally so
### AFD heatmaps are nice
#########################################

library(fields)

u <- 0.05 #bin size
bin.sfs <- read.table(file="../003.output/SFS.derived.form.pub.0.05bin.txt")
bin.sfs <-as.matrix(bin.sfs)
lbs <- log(bin.sfs)
lbs[nrow(lbs),ncol(lbs)] <- NA  ### set this manually, otherwise doesn't work
pdf(file="../003.output/SFS.derived.form.pub.rainbow.0.05bin.pdf", width=6, height=6)
image.plot(lbs,xlab="A. pubescens derived allele frequency", ylab="A. formosa derived allele frequency",legend.lab="number of SNPs (log)")
dev.off()

u <- 0.1 #bin size
bin.sfs <- read.table(file="../003.output/SFS.derived.form.pub.0.1bin.txt")
bin.sfs <-as.matrix(bin.sfs)
lbs <- log(bin.sfs)
lbs[nrow(lbs),ncol(lbs)] <- NA  ### set this manually, otherwise doesn't work
pdf(file="../003.output/SFS.derived.form.pub.rainbow.0.1bin.pdf", width=6, height=6)
image.plot(lbs,xlab="A. pubescens derived allele frequency", ylab="A. formosa derived allele frequency",legend.lab="number of SNPs (log)")
dev.off()


#### could also do SFS of significant derived SNPs
#### should output chr and pos of SNPs that can be polarized





####################################################################################
################## function development below this line ############################
################## function is also updated in source file #########################
####################################################################################
get.derived.af <- function(af.line){
	#print(af.line)
	## parse out allele frequencies
	alleles <- unlist(af.line[5:14])
	alleles <- strsplit(alleles, "/")
	alleles <- do.call(rbind, alleles)
	class(alleles) <- "numeric"
	colnames(alleles) <- unlist(strsplit(as.character(af.line[4]),"/"))
	
	## decide which is derived
	sa <- as.character(af.line[19])
	alleles <- alleles[,colnames(alleles)!=sa]
	## check that semi allele is one of the AQ alleles
	## if not, return NAs
	if(length(alleles)!=10) {out.freqs <- rep(NA,10)
		} else {out.freqs <- alleles}
	return(out.freqs)
}




