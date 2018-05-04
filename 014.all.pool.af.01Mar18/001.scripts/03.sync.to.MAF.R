#### parse filtered sync file to get minor allele counts and minor allele frequences for polymorphic sitesF.helper.fxns.R
#### DLF 05MAR18


InDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input"
InDataFile="no169.pools.poly.sync"
OutDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output"
min.MAC=2
FunFile="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/04.AF.helper.fxns.R"
BamDir="/lustre/scratch/projects/aquilegia/001.pool.bams"
BamList="no169.pool.nodups.bam.list"


setwd(InDataDir)

source(FunFile)


#################################
### step one: parse sync file ###
#################################

all.dat <- read.table(InDataFile,stringsAsFactors=FALSE)

all.dat.freq <- apply(all.dat,1, function(x){parse.line(x,2)})

#test <- apply(all.dat[1:100000,],1,function(x){parse.line(x,2)})

save(all.dat.freq, file=paste(OutDataDir,"all.dat.freq.Rdata", sep="/"))

all.dat.freq <- t(all.dat.freq)


#################################################################################
### step two: only consider sites that are polymorphic after removing low MAC ###
#################################################################################

allele.table <- table(all.dat.freq[,4])
allele.table <- allele.table[order(allele.table)] #about 13 million sites are now non-polymoprhic (i.e. min allele didn't occur frequently enough.  with minMAC=2, this means that these were all singletons.)

poly.freq <- all.dat.freq[nchar(all.dat.freq[,4])>1,]  #22 million to 9 million sites to consider
poly.allele.table <- table(poly.freq[,4])
poly.allele.table <- poly.allele.table[order(poly.allele.table)]
biallele.table <- poly.allele.table[nchar(names(poly.allele.table))==3]  ### 8.3 million are bialleleic (also includes deletions)

biallele.freq <- poly.freq[nchar(poly.freq[,4])==3,]
save(biallele.freq, file=paste(OutDataDir,"biallele.freq.Rdata", sep="/")

#############################################
### step three: AFD between species pools ###
#############################################

# in 05.AF.to.AFD.R





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
		}
	else{	
		t1 <- t.test(alleles[,1]~alleles[,3], data=alleles)
		t1.sum <- c(t1$estimate, t1$p.value)
		}
	names(t1.sum) <- c("f1","p1","pval1")
	
	tt2 <- try(t.test(alleles[,2]~alleles[,3], data=alleles))
	if(length(tt2)==1){
		af2 <- aggregate(alleles[,2]~alleles[,3], data=alleles, mean)
		t2.sum <- c(af2[,2], NA)
		}
	else{
		t2 <- t.test(alleles[,2]~alleles[,3], data=alleles)
		t2.sum <- c(t2$estimate, t2$p.value)
		}
	names(t2.sum) <- c("f2","p2","pval2")
	
	out.dat <- unlist(c(t1.sum, t2.sum))
	return(out.dat)
}


##########################################
### step four: assess with permutation ###
##########################################






##############################
### step five: plot output ###
##############################














##################################### function development below this line ############################################





parse.line <- function(snp.line, min.MAC){
	#print(snp.line)
	
	gts <- unlist(snp.line[4:length(snp.line)])
	gt.m <- sapply(gts, function(x) {strsplit(x, ":")})
	gt.m <- do.call(rbind,gt.m)
	colnames(gt.m) <- c("A","T","C","G","N","D")
	class(gt.m) <- "numeric"
	gt.sums <- apply(gt.m,2,sum)
	## check if min allele occurs at least min.MAC
	gt.check <- gt.sums[gt.sums>=min.MAC]
	if (length(gt.check)>1){
		# get variant bp
		gt.g <- gt.m[,gt.sums>=min.MAC]
		# get AF
		pool.sum <- apply(gt.m,1,sum)
		pool.freq <- apply(gt.g,2,function(x){round(x/pool.sum,digits=3)})
		
		# output variables
		var.bp <- paste(colnames(pool.freq), collapse="/")
		var.freq <- apply(pool.freq,1, function(x){paste(x, collapse="/")})
		out.dat <- c(snp.line[1:3], var.bp, var.freq)
		out.dat <- unlist(out.dat)
	}else {
		gt.g <- gt.check
		pool.sum <- apply(gt.m, 1, sum)
		pool.freq <- gt.m[,colnames(gt.m)==names(gt.check)]/pool.sum  ### won't necessarily add up to 1 since minor alleles are sometimes filtered out if occur less than cutoff
		
		var.bp <- names(gt.g)
		var.freq <- round(pool.freq,3)
		#print(var.freq)
		out.dat <- c(snp.line[1:3], var.bp, var.freq)
		out.dat <- unlist(out.dat)
	}
	return(out.dat)
	
}
