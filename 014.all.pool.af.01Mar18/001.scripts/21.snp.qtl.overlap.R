#### script to look at overlap between sig diff snps and QTL
#### using genome rotation
#### DLF 02 May 18
#### updated 04 May to included fixed SNPs as well


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
qtl <- read.table("./QTL.locations.txt",stringsAsFactors=FALSE,header=TRUE)

###3. SNP allele frequencies
load("./after.perm.Rdata")
#both.mean.af

load("./perm.results.Rdata")
sig.snps  <- perm.results[perm.results[,3]==0 & perm.results[,4]==0 & perm.results[,5]==0 & perm.results[,6]==0,]
sig.snps <- sig.snps[is.na(sig.snps$p.extreme)==FALSE,]


#############################################
#### change variables to absolute base pair index
#############################################
### qtl interval positions
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

### qtl peak positions
### taking midpoint of bin +/-250kb
### it looks like that is the overall min bin size???
qtl.peak <- cbind(qtl.n[,1]-250000, qtl.n[,1]+250000


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
ss.pos <- ss.n[,2]


######################################################
#### rotation enrichment for different pvalue cutoffs and determinations
#######################################################
## proportion of sig snps overlapping a qtl
## ss.pos is position of sig. differentiated SNPs (however you define) - a vector
## qtl.pos is start and stop position of QTL peaks - a numerical matrix
#ss.pos <- ss.n[,2]
#qtl.pos <- qtl.n[,2:3]
## n.rotations is the number of rotations to do
## max.bp is the length in bp of the genome

### trying SNPS with t-test p-value less than bonferroni
### also include fixed differences as defined previously
tp <- af.n[af.n$fp.tt<=(0.05/nrow(af.n)) | is.na(af.n$fp.cons)==FALSE,]
tp <- tp[is.na(tp[,1])==FALSE,]  ###100 SNPs
tp.pos <- as.numeric(tp[,2])
bonf.num <- qtl.snp.rotation(ss.pos=tp.pos, qtl.pos=qtl.n[,2:3], max.bp=len.max, n.rotation=1000) ##308
bonf.num.peak <- qtl.snp.rotation(ss.pos=tp.pos, qtl.pos=qtl.peak, max.bp=len.max, n.rotation=1000) ##1000


### using diverged snps according to permutations
perm.num <- qtl.snp.rotation(ss.pos=ss.n[,2], qtl.pos=qtl.n[,2:3], max.bp=len.max, n.rotation=1000) #121
perm.num.peak <- qtl.snp.rotation(ss.pos=ss.n[,2], qtl.pos=qtl.peak, max.bp=len.max, n.rotation=1000) #121


## t-test p-value <=0.000001 (-log10 ==6)
### also include fixed differences as defined previously
tp <- af.n[af.n$fp.tt<=0.000001 | is.na(af.n$fp.cons)==FALSE,]
tp <- tp[is.na(tp[,1])==FALSE,]  ###4109 SNPs
tp.pos <- as.numeric(tp[,2])
pval6.num <- qtl.snp.rotation(ss.pos=tp.pos, qtl.pos=qtl.n[,2:3], max.bp=len.max, n.rotation=1000) #106
pval6.num.peak <- qtl.snp.rotation(ss.pos=tp.pos, qtl.pos=qtl.peak, max.bp=len.max, n.rotation=1000) #939


tp <- af.n[af.n$fp.tt<=0.00001 | is.na(af.n$fp.cons)==FALSE,]
tp <- tp[is.na(tp[,1])==FALSE,] #13728
tp.pos <- as.numeric(tp[,2])
pval5.num <- qtl.snp.rotation(ss.pos=tp.pos, qtl.pos=qtl.n[,2:3], max.bp=len.max, n.rotation=1000) #128
pval5.num.peak <- qtl.snp.rotation(ss.pos=tp.pos, qtl.pos=qtl.peak, max.bp=len.max, n.rotation=1000) #902


save.image("../004.rotation.output/21.snp.qtl.overlap.Rdata")




############ testing below this line ######################################
qtl.snp.rotation <- function(ss.pos, qtl.pos, max.bp, n.rotations){
        ### get observed values
        obs.prop <- snp.qtl.prop(ss.pos=ss.pos, qtl.pos=qtl.pos)

        ### get rotated values
        rot.bp <- sample(1:max.bp, n.rotations, replace=FALSE)
        #print(rot.bp)
        rot.prop <- sapply(rot.bp, function(bp){
        #       print(bp)
                pos.r <- genome.rotate(pos=ss.pos, bp.slide=bp, max.bp=len.max)
                prop.r <- snp.qtl.prop(ss.pos=pos.r, qtl.pos=qtl.pos)
                return(prop.r)
                })
        rot.extreme <- rot.prop[rot.prop>=obs.prop]
        n.extreme <- length(rot.extreme)
        return(n.extreme)
        }
