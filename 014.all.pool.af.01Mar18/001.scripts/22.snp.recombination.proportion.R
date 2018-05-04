#### look at proportion of snps in each recombination rate bin
#### rotate genome to look for significance
#### DLF 02MaY18


setwd("/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input")
source("/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/20.genome.rotation.functions.R")


##########################
#### get input variables
###########################

###1. length of chromosomes
len <- read.table("/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta.fai",stringsAsFactors=FALSE, nrows=7)
len.cs <- cumsum(len[,2])
len.cs <- c(0,len.cs)
len.max <- max(len.cs)

###2. positions of QTLs
### don't need qtl at this point, but might later, so just keep.
qtl <- read.table("./QTL.locations.txt",stringsAsFactors=FALSE,header=TRUE)

###3. recombination rate in bins
recomb <- read.table("/lustre/scratch/projects/aquilegia/007.evol.analysis/008.recombination/recomb.pi.data.out.txt",stringsAsFactors=FALSE, header=TRUE)
recomb$midpt <- recomb$starts+500000.5
### step one: set H,M,L borders
quantile(recomb$p.cm, na.rm=TRUE)
r.nonzero <- recomb[recomb$p.cm>0,]
quantile(r.nonzero$p.cm, na.rm=TRUE)
### so... let's say low is 0, medium is >0 to 1.575 (which is 50th percentile of nonzero values), high is > 1.575
rec.windows <- recomb[,c(1,2,23,13)]
rec.bin <- apply(rec.windows,1, function(x){
        up.cm <- as.numeric(x[3])
        bin <- "e"
        if(is.na(up.cm)) {bin <- NA
        } else if (up.cm==0) {bin <- "l"
        } else if (up.cm>1.575) {bin <- "h"
        } else {bin <- "m"}
        return(bin)
})
rec.windows$rec.bin <- rec.bin
rec.windows$rw.id <- c(1:nrow(rec.windows))

###4. SNP allele frequencies
load("./after.perm.Rdata")
both.mean.af <- cbind(both.mean.af, fp.tt)
load("./perm.results.Rdata")
sig.snps  <- perm.results[perm.results[,3]==0 & perm.results[,4]==0 & perm.results[,5]==0 & perm.results[,6]==0,]
sig.snps <- sig.snps[is.na(sig.snps$p.extreme)==FALSE,]


#############################################
#### change variables to absolute base pair index
#############################################
### qtl positions
qtl.s <- qtl[,c(3,5,12,13)]
qtl.s <- split(qtl.s, qtl.s$chr)
for (up in 1:7){
	uq <- qtl.s[[up]]
	ul <- len.cs[up]
	uq <- uq[,2:4]*1000000
	uq <- apply(uq,2, function(x){x+ul})
	qtl.s[[up]]<- uq
}
qtl.n <- do.call(rbind, qtl.s)

### recomb windows
rw.s <- split(rec.windows, rec.windows$chr)
for(up in 1:7){
	ur <- rw.s[[up]]
	ul <- len.cs[up]
	ur[,1] <- ur[,1]+ul
	ur[,2] <- ur[,2]+ul
	rw.s[[up]] <- ur
}
rw.n <- do.call(rbind, rw.s)

### SNP allele frequencies and pvals
af.s <- split(both.mean.af,both.mean.af[,1])
for(up in 1:7){
	uaf <- af.s[[up]]
	ul <- len.cs[up]
	uaf[,2] <- uaf[,2]+ul
	af.s[[up]] <- uaf
}
af.n <- do.call(rbind, af.s)
rownames(af.n)<- c(1:nrow(af.n))

ss.s <- split(sig.snps, sig.snps$chr)
for(up in 1:7){
	us <- ss.s[[up]]
	ul <- len.cs[up]
	us[,2] <- us[,2]+ul
	ss.s[[up]] <- us
}
ss.n <- do.call(rbind, ss.s)
ss.n <- ss.n[,1:2]
ss.pos <- as.numeric(ss.n[,2])

############################################
#### do rotation
################################################

### sig divergence snps by permutation

t.srn <- snp.recomb.num(ss.pos=ss.pos, rw.n=rw.n)

rot.bp <- sample(1:(len.max), 1000, replace=FALSE)
rot.srn <- sapply(rot.bp, function(bp){
	print(bp)
	pos.r <- genome.rotate(pos=ss.pos, bp.slide=bp, max.bp=len.max)
	num.srn <- snp.recomb.num(ss.pos=pos.r, rw.n=rw.n)
	return(num.srn)
})

save(rot.srn, file="../004.rotation.output/rot.srn.Rdata")

### what to do with this now???
### then write function to do this
### and test other pvalue cutoffs


