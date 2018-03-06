#### parse filtered sync file to get minor allele counts and minor allele frequences for polymorphic sites
#### DLF 05MAR18


InDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input"
InDataFile="no169.pools.poly.sync"
OutDataDir="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output"
min.MAC=2
FunFile="/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/sync.to.MAF.helper.fxns.R"



setwd(InDataDir)

load(FunFile)


#################################
### step one: parse sync file ###
#################################

all.dat <- read.table(InDataFile,stringsAsFactors=FALSE)

all.dat.freq <- apply(all.dat,1, function(x){parse.line(x,2)})

#test <- apply(all.dat[1:100000,],1,function(x){parse.line(x,2)})

save(all.dat.freq, file=paste(OutDataDir,"all.dat.freq.Rdata", sep="/"))


#################################################################################
### step two: only consider sites that are polymorphic after removing low MAC ###
#################################################################################






#############################################
### step three: AFD between species pools ###
#############################################



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