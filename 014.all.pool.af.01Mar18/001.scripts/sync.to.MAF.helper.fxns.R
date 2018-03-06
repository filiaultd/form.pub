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
