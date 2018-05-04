### parse allele frequencies of bialleic snps
### based off similar analysis of fst by gene
### DLF 17April18

input.dir <- "/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input/"
output.dir <- "/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output/"

gs.file <- "no169.pools.poly.sync.goodsites.out"
bi.file <- "no169.pools.poly.sync.biallelic.out"


setwd(input.dir)

#bi.snps <- read.table(bi.file, stringsAsFactors=FALSE)
#chrs <- unique(bi.snps[,1])[1:7]
#bi.snps <- bi.snps[bi.snps[,1]%in%chrs,]

#pool.index <- read.table("pool.order.index.txt", header=TRUE,stringsAsFactors=FALSE) 

#load(file="biallelic.allele.freq.bi.af.Rdata")  ### this is allele frequencies

##get sig diffs, plot across genome, permute as others

#pub.mean.af <- apply(bi.af[,which(pool.index$species=="pubescens")], 1, function(x){mean(x, na.rm=TRUE)})
#form.mean.af <- apply(bi.af[,which(pool.index$species=="formosa")], 1, function(x){mean(x, na.rm=TRUE)})

#mean.fp.diff <- form.mean.af-pub.mean.af
#both.mean.af <- cbind(bi.snps[,1:2],pub.mean.af, form.mean.af, mean.fp.diff)




## can reuse this function although the variable names aren't very fitting now...  ;)

fst.kw <- function(fst.a, fst.b){
	pval.kw <- rep(NA,nrow(fst.a))
	for(up in 1:nrow(fst.a)){
		a.dat <- fst.a[up,]
		a.dat <- a.dat[is.na(a.dat)==FALSE]
		b.dat <- fst.b[up,]
		b.dat <- b.dat[is.na(b.dat)==FALSE]
		kwt <- try(kruskal.test(list(a.dat, b.dat), silent=TRUE))
		if (length(kwt)>1){pval.kw[up]<- kwt$p.value}
	}
	return(pval.kw)
}


#fp.kw <- fst.kw(bi.af[,which(pool.index$species=="formosa")],bi.af[,which(pool.index$species=="pubescens")])

#both.mean.af <- cbind(both.mean.af, fp.kw)

#save.image()
load(".RData")


#pdf("../003.output/mean.allele.frequencies.pdf", width=5, height=5)
#hist(abs(both.mean.af$mean.fp.diff), xlab="mean allele freq difference", main="Histogram of AFD between form and pub")
#dev.off()

### nothing significant with KW test!
### could try t-test?

############################################
##### plot AFDs across genome
############################################

### first set variables for genome
gen.pos <- read.table("/lustre/scratch/projects/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta.fai",stringsAsFactors=FALSE)[1:7,]
max.pos <- gen.pos[,2]
cs <- cumsum(max.pos)
bp.add <- c(0, cs)[1:7]
gl <- sum(max.pos)
c.col <- c("darkblue","dodgerblue1","darkblue", "dodgerblue1", "darkblue","dodgerblue1","darkblue", "dodgerblue1")
cp.start <- bp.add[1:7]+(0.5*max.pos)

### var.label is the name of the variable for labeling plot y axis
### win.var is a data frame of chr pos value pvalue(-log10)
### plot.name is main text for plot
### pval.cutoff is obvious
### cex.plot (change if needed for size reasona)

genome.plot.pval <- function(var.label, win.var, plot.name, cex.plot=1, pval.cutoff){
	up.dat <- split(win.var, win.var$chr)
	max.y <- max(win.var[,3], na.rm=TRUE)
    min.y <- min(win.var[,3], na.rm=TRUE)
    plot(0,0, xlim=c(0,gl),ylim=c(min.y,max.y),type="n", xaxt="n",ylab="",xlab="")
    for(up.c in 1:7){
        c.dat <- up.dat[[up.c]][,3]
        up.pos <- up.dat[[up.c]][,2]
        c.bp.add <- bp.add[up.c]
        c.pos <- up.pos + c.bp.add
        points(c.pos, c.dat,type="p", col=c.col[up.c], pch=19)
        ### highlight sig
        sig.dat <-up.dat[[up.c]]
        sig.dat <- sig.dat[sig.dat[,4]>pval.cutoff,]
        sig.dat <- sig.dat[is.na(sig.dat[,1])==FALSE,]
        s.pos <- sig.dat[,2] + c.bp.add
        s.point <- sig.dat[,3]
        points(s.pos, s.point, type="p", col="purple", pch=19)
        	}
        abline(v=c(0,cs), col="grey50",lty=3)
        mtext(plot.name, side=3, line=-1.5, cex=cex.plot, adj=0.975, font=4)
        mtext(1:7, 1,at=cp.start, line=0.5, cex=cex.plot)
        mtext("Chromosome", 1, line=2, cex=cex.plot)
        mtext(var.label,2,las=3,cex=cex.plot, line=2.2)
}


#pdf("../003.output/fst.differences.with.pvals.pdf", width=20, height=6)
#win.var <- both.mean.af[,c(1,2,5,6)]
#win.var <- win.var[abs(win.var$mean.fp.diff)>0.6,]
#win.var[,4] <- -log10(win.var[,4])
#win.var[,3] <- abs(win.var[,3])
#colnames(win.var)[1:2]<-c("chr","pos")
#pval.cutoff=-log10(0.05/nrow(win.var))
#genome.plot.pval(var.label="absolute AFD between species", win.var=win.var, plot.name="", pval.cutoff=pval.cutoff)	
#dev.off()

##########################################################
#### what do AFDs look like if you do pool ID permutation?
############################################################


#### permute one pool per species

#### permute one sample per species
po.s <- split(pool.index, pool.index$species)

one.switch <- expand.grid(po.s[[1]][,3], po.s[[2]][,3])
colnames(one.switch) <- c("form","pub")
afd.perm.out <- matrix(NA, nrow=nrow(bi.af), ncol=nrow(one.switch))

for(up in 1:nrow(one.switch)){
	print(up)
	po.p <- pool.index
	po.p[one.switch[up,1],5] <- "pubescens"
	po.p[one.switch[up,2],5] <- "formosa"
	pp.mean.af <- apply(bi.af[,which(po.p$species=="pubescens")], 1, function(x){mean(x, na.rm=TRUE)})
	pf.mean.af <- apply(bi.af[,which(po.p$species=="formosa")], 1, function(x){mean(x, na.rm=TRUE)})
	mean.perm.diff <- pf.mean.af-pp.mean.af
	afd.perm.out[,up] <- mean.perm.diff
}

### permute two samples per species
two.switch <- expand.grid(po.s[[1]][,3], po.s[[1]][,3], po.s[[2]][,3], po.s[[2]][,3])
two.switch <- two.switch[two.switch[,1]<two.switch[,2] & two.switch[,3]<two.switch[,4],]
afd.perm2.out <- matrix(NA, nrow=nrow(bi.af), ncol=nrow(two.switch))
for(up in 1:nrow(two.switch)){
	print(up)
	po.p <- pool.index
	po.p[two.switch[up,1],5] <- "pubescens"
	po.p[two.switch[up,2],5] <- "pubescens"
	po.p[two.switch[up,3],5] <- "formosa"
	po.p[two.switch[up,4],5] <- "formosa"
	pp.mean.af <- apply(bi.af[,which(po.p$species=="pubescens")], 1, function(x){mean(x, na.rm=TRUE)})
	pf.mean.af <- apply(bi.af[,which(po.p$species=="formosa")], 1, function(x){mean(x, na.rm=TRUE)})
	mean.perm.diff <- pf.mean.af-pp.mean.af
	afd.perm2.out[,up] <- mean.perm.diff	
}

save.image("redone.perms.Rdata")

######################################################
### how often do you get AFD with the same profile as with the "real" split?
### how to judge this?
################################################

#perm1.prop <- matrix(NA, nrow=nrow(both.mean.af), ncol=2)
#for (up in 1:nrow(both.mean.af)){
#	o.dat <- both.mean.af[up,5]
#	p.dat <- afd.perm.out[up,]
#	p.dat <- p.dat[is.na(p.dat)==FALSE]
#	n.perm <- length(p.dat)
#	### how often is permuted allele frequency diff < or equal to observed
#	if(is.na(o.dat)==TRUE){
#		p.extreme <- NA
#		p.eq <- NA
#	} else if(o.dat <0){
#		p.extreme <- length(p.dat[p.dat<o.dat])/n.perm
#		p.eq <- length(p.dat[p.dat==o.dat])/n.perm
#	} else {
#		p.extreme <- length(p.dat[p.dat>o.dat])/n.perm
#		p.eq <- length(p.dat[p.dat==o.dat])/n.perm
#	}
#	n.out <- c(p.extreme, p.eq)
#	perm1.prop[up,] <- n.out
#}
#colnames(perm1.prop) <-c("p.extreme", "p.eq")

#perm2.prop <- matrix(NA, nrow=nrow(both.mean.af), ncol=2)
#for (up in 1:nrow(both.mean.af)){
#	o.dat <- both.mean.af[up,5]
#	p.dat <- afd.perm2.out[up,]
#	p.dat <- p.dat[is.na(p.dat)==FALSE]
#	n.perm <- length(p.dat)
#	### how often is permuted allele frequency diff < or equal to observed
#	if(is.na(o.dat)==TRUE){
#		p.extreme <- NA
#		p.eq <- NA
#	} else if(o.dat <0){
#		p.extreme <- length(p.dat[p.dat<o.dat])/n.perm
#		p.eq <- length(p.dat[p.dat==o.dat])/n.perm
#	} else {
#		p.extreme <- length(p.dat[p.dat>o.dat])/n.perm
#		p.eq <- length(p.dat[p.dat==o.dat])/n.perm
#	}
#	n.out <- c(p.extreme, p.eq)
#	perm2.prop[up,] <- n.out
#}
#colnames(perm2.prop) <-c("p.extreme2", "p.eq2")

#perm.results <- cbind(bi.snps[,1:2], perm1.prop, perm2.prop)
#colnames(perm.results)[1:2] <- c("chr", "pos")
#save(perm.results, file="perm.results.Rdata")

############################################################################
### take only SNPs where no permutation is more extreme or equal to observed value - most conservative I can be!
###################################################
#sig.snps  <- perm.results[perm.results[,3]==0 & perm.results[,4]==0 & perm.results[,5]==0 & perm.results[,6]==0,]
####  449546 SNPs!  so many!
#### what do these AFDs look like?

#colnames(both.mean.af)[1:2] <- c("chr","pos")

#sig.af <- merge(sig.snps, both.mean.af, all.x=TRUE)

### compare ecdf of all snps vs "sig" snps
#all.ecdf <- ecdf(abs(both.mean.af$mean.fp.diff))
#sig.ecdf <- ecdf(abs(sig.af$mean.fp.diff))

#pdf("../003.output/ecdf.all.vs.sig.snps.pdf", width=5, height=5)
#plot(all.ecdf, do.points=FALSE, verticals=TRUE, col="green", xlab="allele frequency difference", ylab="cumulative proportion")
#plot(sig.ecdf, do.points=FALSE, verticals=TRUE, add=FALSE, col="blue")
#legend(x="topright", col=c("green", "blue"),legend=c("all snps","sig snps"))
#dev.off()

### can aslo plot regular histograms
#pdf("../003.output/histograms.AFD.all.and.sig.snps.pdf")
#hist(both.mean.af$mean.fp.diff, main="all snps", xlab="allele frequency difference")
#hist(sig.af$mean.fp.diff, main="sig snps", xlab="allele frequency difference")
#dev.off()

#save(sig.af, file="../003.output/sig.snps.Rdata") ### so can use this in overlap script

### som some of the intermediate files are missing here due to accidental loss of .Rdata file.  Can redo, but I have the output files anyway, so....

################################################
### another thing to try is to just do a regular t-test 
### rather than KW
#################################################
#load("after.perm.Rdata")

afd.ttest <- function(afd.a, afd.b){
	pval.tt <- rep(NA,nrow(afd.a))
	for(up in 1:nrow(afd.a)){
		a.dat <- afd.a[up,]
		a.dat <- a.dat[is.na(a.dat)==FALSE]
		b.dat <- afd.b[up,]
		b.dat <- b.dat[is.na(b.dat)==FALSE]
		tt <- try(t.test(a.dat, b.dat), silent=TRUE)
		if (length(tt)>1){pval.tt[up]<- tt$p.value}
	}
	return(pval.tt)
}


fp.tt <- afd.ttest(bi.af[,which(pool.index$species=="formosa")],bi.af[,which(pool.index$species=="pubescens")])

# the problem here, though, is that sometimes the data are consistant (i.e.  all 1 in one species and all 0 in the other)
# we should also get an idea about those guys, too.

afd.cons <- function(afd.a, afd.b){
	fix.afd <- rep(NA, nrow(afd.a))
	for (up in 1:nrow(afd.a)){
		a.dat <- afd.a[up,]
		au <- length(unique(a.dat))
		b.dat <- afd.b[up,]
		bu <- length(unique(b.dat))
		if(au==1 & bu==1){fix.afd[up] <- paste(au, bu, sep="/")}
	}
	return(fix.afd)
} 

fp.cons <- afd.cons(bi.af[,which(pool.index$species=="formosa")],bi.af[,which(pool.index$species=="pubescens")])
## 35 SNPs aboslutely split
### let's just make sure these aren't missing a lot of data
cons.dat <- bi.af[rownames(both.mean.af[is.na(both.mean.af$fp.cons)==FALSE,]),]
cd.tab <- table(is.na(cons.dat))  ### data present for each one.
### should also note at AFD is 1 for all of these (except one where it is 0.5) ( Chr_06 19568838 )

both.mean.af <- cbind(both.mean.af,fp.tt)
both.mean.af <- cbind(both.mean.af, fp.cons)
save(both.mean.af, file="both.mean.af.Rdata")


#save.image("after.perm.Rdata")
load("after.perm.Rdata")

bp.tmp <- 0.05/length(fp.tt)
sig.ttest <- both.mean.af[both.mean.af$fp.tt<=bp.tmp,]
sig.ttest <- sig.ttest[is.na(sig.ttest[,1])==FALSE,]  ### 65 snps pass bonferroni on all chromosomes except 4...
write.table(sig.ttest, file="../003.output/ttest.sig.snps.txt", quote=FALSE)


######################################################
### genome rotation for overlap will be
### in another script
#####################################################


#######################################################
#### do some plotting of significant AFDs across the genome
##########################################################
###########################
### work locally 
### this means reloading output files
############################


setwd("/Volumes/aquilegia/004.pooled.sequencing/014.all.pool.af.01Mar18/002.input/")

load("both.mean.af.Rdata")
load("../003.output/sig.snps.Rdata")
colnames(both.mean.af)[1:2] <- c("chr","pos")

bma.s <- split(both.mean.af, both.mean.af$chr)
ss.s <- split(sig.af, sig.af$chr)

for(up in 1:7){
	up.b <- bma.s[[up]]
	up.s <- ss.s[[up]]
	up.out <- merge(up.s, up.b, by="pos")
	ss.s[[up]] <- up.out
}

sig.dat <- do.call(rbind, ss.s)


##### two plots:
##### 1.  Plot AFD with p-values marked across genome

### var.label is the name of the variable for labeling plot y axis
### win.var is a data frame of chr pos value pvalue(-log10) p.cons(which is whether or not data is consistant within both pools)
### plot.name is main text for plot
### pval.cutoff is obvious
### cex.plot (change if needed for size reasons)
#### recomb is my recombination dataset

genome.plot.pval <- function(var.label, win.var, plot.name, cex.plot=1, pval.cutoff, recomb){
	par(mar=c(4,5,1,5))
	up.dat <- split(win.var, win.var$chr)
	max.y <- max(win.var[,3], na.rm=TRUE)
    min.y <- min(win.var[,3], na.rm=TRUE)
    plot(0,0, xlim=c(0,gl),ylim=c(min.y,max.y),type="n", xaxt="n",ylab="",xlab="")
    for(up.c in 1:7){
        c.dat <- up.dat[[up.c]][,3]
        up.pos <- up.dat[[up.c]][,2]
        c.bp.add <- bp.add[up.c]
        c.pos <- up.pos + c.bp.add
        points(c.pos, c.dat,type="p", col=c.col[up.c], pch=19)
        ### highlight sig
        sig.dat <-up.dat[[up.c]]
        sig.dat <- sig.dat[sig.dat[,4]>pval.cutoff,]
        sig.dat <- sig.dat[is.na(sig.dat[,1])==FALSE,]
        s.pos <- sig.dat[,2] + c.bp.add
        s.point <- sig.dat[,3]
        points(s.pos, s.point, type="p", col="purple", pch=19)
        ### also highlight invariable
        in.dat <- up.dat[[up.c]]
        in.dat <- in.dat[is.na(in.dat$fp.cons)==FALSE,]
        in.pos <- in.dat[,2] + c.bp.add
        in.point <- in.dat[,3]
        points(in.pos, in.point, type="p", col="red", pch=19)
        }
      ## plot recombination on second axis
      rec <- recomb
      rec.s <- split(rec, rec$chr)
      par(new=T)
      plot(1, 1, type="n", xlim=c(0,gl), ylim=c(0,max(rec$p.cm, na.rm=TRUE)), xaxt="n", yaxt="n",xlab="", ylab="")
      for(up in 1:7){
      	c.bp.add <- bp.add[up]
      	up.r <- rec.s[[up]]
     	r.pos <- up.r$midpt + c.bp.add
     	r.point <- up.r$p.cm
     	points(r.pos, r.point, type="l", col="darkblue", lwd=3)
      }
     axis(side=4, cex=cex.plot)
        ### annotate plot
        abline(v=c(0,cs), col="grey50",lty=3)
        mtext(plot.name, side=3, line=-1.5, cex=cex.plot, adj=0.975, font=4)
        mtext(1:7, 1,at=cp.start, line=0.5, cex=cex.plot)
        mtext("Chromosome", 1, line=2, cex=cex.plot)
        mtext(var.label,2,las=3,cex=cex.plot, line=2.5)
        mtext("cM/Mbp",4,las=3,cex=cex.plot, line=2.5)
        #legend(x="bottomright",col=c("purple", "red"), legend=c(paste("sig at -log10 pval=", pval.cutoff, sep=""), "consistant within species"), pch=19, bg="white", cex=1.2)
}



recomb <- read.table("/Volumes/aquilegia/002.species.mapping/007.evol.analysis/008.recombination/recomb.pi.data.out.txt",stringsAsFactors=FALSE, header=TRUE)
recomb$midpt <- recomb$starts+500000.5


win.var <- sig.dat[,c(2,1,9,16,17)]
win.var[,3] <- abs(win.var[,3])
win.var[,4] <- -log10(win.var[,4])

gen.pos <- read.table("/Volumes/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta.fai",stringsAsFactors=FALSE)[1:7,]
max.pos <- gen.pos[,2]
cs <- cumsum(max.pos)
bp.add <- c(0, cs)[1:7]
gl <- sum(max.pos)
#c.col <- c("darkblue","dodgerblue1","darkblue", "dodgerblue1", "darkblue","dodgerblue1","darkblue", "dodgerblue1")
c.col <- c("grey75","grey55","grey75", "grey55", "grey75","grey55","grey75")
cp.start <- bp.add[1:7]+(0.5*max.pos)


pdf("../003.output/allele.freq.diff.with.recomb.pdf", width=20, height=8)
#quartz(width=20, height=8)
genome.plot.pval(var.label="allele frequency difference", win.var=win.var, plot.name="",cex.plot=1.2, pval.cutoff=7, recomb=recomb)
dev.off()


##### 2.  Plot joint SFS of significant SNPs.
##### do this once SNPs have been properly polarized.



