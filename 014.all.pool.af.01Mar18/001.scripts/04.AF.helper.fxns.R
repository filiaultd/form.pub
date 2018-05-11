### helper functions for R scripts to look at allele frequency differences in pools
### DLF 06March18



#### parse a line of a sync file(snp.line), filter for minimim allele cound (min.MAC), output chr, pos, ref, alleles, and allele frequencies by site for all alleles

parse.line <- function(snp.line, min.MAC){
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
		}
	else {
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



###### from a line of biallelic allele frequencies, calculate AFD and p-value for both alleles
###### comparing form and pub pools
######

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


###### same as get.af, but polarize SNP to Semi
##### returns frequency and p-val of derived variants
##### used in 25.polarization.R - see this script for input file format

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






