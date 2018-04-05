### parse FST by gene data
### DLF 21 March 18

setwd("/Volumes/aquilegia/004.pooled.sequencing/016.popool2.all.pops.16March18/003.output/")
data.file <- 'no169.genes.15MAR18.fst'
gtf.file <- ("/Volumes/aquilegia/001.common.reference.files/018.Acoerulea_v3.1_annotation_20150930/Acoerulea_396_v3.1_genes_only.gtf")
outdir <- "/Volumes/aquilegia/004.pooled.sequencing/016.popool2.all.pops.16March18/004.parsed.data/"

###############################################
### step one: enter data and format variables
################################################
fst <- read.table(file=data.file, stringsAsFactors=FALSE)
p.comb <- unlist(fst[1,6:50])
p.comb <- unlist(sapply(p.comb, function(x){strsplit(x, "=")[[1]][1]}))
colnames(fst) <- c("gene","window","snp.no","prop","something",p.comb)
fst <- fst[,c(1,6:50)]

### parse to get chromosome and positions

gtf <- read.table(gtf.file, stringsAsFactors=FALSE)
gtf <- gtf[,c(1,4,5,10)]
colnames(gtf) <- c("chrom", "start", "stop", "gene")
fst <- merge(gtf,fst, by="gene",all=TRUE)
fst <- fst[grep("Chr", fst$chrom),]  #29550 primary gene models

rm.comb <- function(x){
	x <- unlist(x)
	out.fst <- sapply(x, function(up){
		up.fst <- strsplit(up,"=")[[1]][2]
		up.fst <- as.numeric(up.fst)
		return(up.fst)
	})
}

tmp <- apply(fst[5:ncol(fst)], 2, rm.comb)
fst[5:ncol(fst)] <- tmp

## need index of pool identification

p.order <- read.table("/Volumes/aquilegia/004.pooled.sequencing/002.bwa.mapping.pooled/006.24Mar2015.mapping/001.pool.bams/no169.pool.nodups.bam.list",stringsAsFactors=FALSE)
p.order$order <- c(1:10)
p.order$pool.number <- sapply(p.order[,1],function(x){strsplit(x,split=".",fixed=TRUE)[[1]][1]})
p.ind <- read.csv("/Volumes/aquilegia/004.pooled.sequencing/015.popool.all.pops.16March18/002.input/pool.index.csv",stringsAsFactors=FALSE)
p.order <- merge(p.order, p.ind)
write.table(p.order, "pool.order.index.txt")

## which combos are between and within species

po.s <- split(p.order, p.order$species)
bet.a <- expand.grid(po.s[[1]][,3], po.s[[2]][,3])
bet.b <- expand.grid(po.s[[2]][,3], po.s[[1]][,3])
b.combos <- rbind(bet.a, bet.b)
b.combos <- paste(b.combos[,1],b.combos[,2],sep=":")
f.a <- combn(po.s[[1]][,3],2)
f.b <- f.a[c(2,1),]
f.combos <- cbind(f.a,f.b)
f.combos <- paste(f.combos[1,],f.combos[2,],sep=":")
p.a <- combn(po.s[[2]][,3],2)
p.b <- p.a[c(2,1),]
p.combos <- cbind(p.a,p.b)
p.combos <- paste(p.combos[1,],p.combos[2,],sep=":")


setwd(outdir)

#######################################################
### step two: look at FST between formosa and pubescens pools
#########################################################

b.fst <- fst[,colnames(fst)%in%b.combos]
w.fst <- fst[,colnames(fst)%in%b.combos==FALSE]
w.fst <- w.fst[,-c(1:4)]
p.fst <- fst[,colnames(fst)%in%p.combos]
f.fst <- fst[,colnames(fst)%in%f.combos]

b.mean <- apply(b.fst,1, mean)
w.mean <- apply(w.fst,1, mean)
p.mean <- apply(p.fst,1, mean)
f.mean <- apply(f.fst, 1,mean)

all.fst <- cbind(b.mean, w.mean, p.mean, f.mean)
colnames(all.fst) <- c("between","within.both", "pubescens", "formosa")
pdf("boxplot.mean.fst.pdf", width=6, height=6)
boxplot(all.fst, ylab="mean FST")
dev.off()


##### dot plots of means

library(ggplot2)
theme_set(theme_bw())
### between vs within
up.plot.b <- qplot(b.mean,w.mean,colour=I("navy"),alpha=I(0.2)) + xlab("mean between-species Fst") + ylab("mean within-species Fst") + xlim(0, max(c(max(b.mean, na.rm=TRUE), max(w.mean, na.rm=TRUE)))) + ylim(0, max(c(max(b.mean, na.rm=TRUE), max(w.mean, na.rm=TRUE))))
scatter <- up.plot.b + geom_density2d(colour=I("orchid"))+ geom_rug(col=I("seagreen4"),alpha=I(0.02))+theme_bw(base_size = 21) 
ggsave("mean.Fst.between.within.pdf",family="Helvetica")
### form vs pub
up.plot.b <- qplot(f.mean,p.mean,colour=I("navy"),alpha=I(0.2)) + xlab("mean formosa Fst") + ylab("mean pubescens Fst") + xlim(0, max(c(max(f.mean, na.rm=TRUE), max(p.mean, na.rm=TRUE)))) + ylim(0, max(c(max(f.mean, na.rm=TRUE), max(p.mean, na.rm=TRUE))))
scatter <- up.plot.b + geom_density2d(colour=I("orchid"))+ geom_rug(col=I("seagreen4"),alpha=I(0.02))+theme_bw(base_size = 21) 
ggsave("mean.Fst.form.pub.pdf",family="Helvetica")
### between vs form
up.plot.b <- qplot(b.mean,f.mean,colour=I("navy"),alpha=I(0.2)) + xlab("mean between-species Fst") + ylab("mean within-formosa Fst") + xlim(0, max(c(max(b.mean, na.rm=TRUE), max(f.mean, na.rm=TRUE)))) + ylim(0, max(c(max(b.mean, na.rm=TRUE), max(f.mean, na.rm=TRUE))))
scatter <- up.plot.b + geom_density2d(colour=I("orchid"))+ geom_rug(col=I("seagreen4"),alpha=I(0.02))+theme_bw(base_size = 21) 
ggsave("mean.Fst.between.within.formosa.pdf",family="Helvetica")
### between vs pub
up.plot.b <- qplot(b.mean,p.mean,colour=I("navy"),alpha=I(0.2)) + xlab("mean between-species Fst") + ylab("mean within-pubescens Fst") + xlim(0, max(c(max(b.mean, na.rm=TRUE), max(p.mean, na.rm=TRUE)))) + ylim(0, max(c(max(b.mean, na.rm=TRUE), max(p.mean, na.rm=TRUE))))
scatter <- up.plot.b + geom_density2d(colour=I("orchid"))+ geom_rug(col=I("seagreen4"),alpha=I(0.02))+theme_bw(base_size = 21) 
ggsave("mean.Fst.between.within.pubescens.pdf",family="Helvetica")


##### test for significant differences

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

bw.kw <- fst.kw(b.fst,w.fst)
fp.kw <- fst.kw(p.fst,f.fst)
bf.kw <- fst.kw(b.fst,f.fst)
bp.kw <- fst.kw(b.fst, p.fst)

all.kw <- cbind(fst[,1:4],bw.kw, fp.kw,bf.kw,bp.kw)

pdf("kw.pval.histograms.pdf", width=8, height=8)
par(mfcol=c(2,2))
mains <- c(rep(NA,4),"between_within","form_pub","between_form","between_pub")
for(up in 5:8){
	hist(-log10(all.kw[,up]), main=mains[up], xlab="-log10 pval")
}
dev.off()
### so admittedly there are power issues here, but the between distributions are very different from the form_pub
### will permute to make power issues less of an issue

### bonferroni cutoff for 29550 genes is 5.7715


##### look at FST differences

bw.diff <- b.mean-w.mean
fp.diff <- f.mean-p.mean
bf.diff <- b.mean-f.mean
bp.diff <- b.mean-p.mean

all.diff <- cbind(bw.diff, fp.diff, bf.diff, bp.diff)

##### plot these across genome
##### make fxn to plot from chr pos value so this is modular
### first set variables for genome
gen.pos <- read.table("/Volumes/aquilegia/001.common.reference.files/017.Aquilegia_coerulea_20150331/20150331/sequences/Aquilegia_coerulea.main_genome.scaffolds.fasta.fai",stringsAsFactors=FALSE)[1:7,]
max.pos <- gen.pos[,2]
cs <- cumsum(max.pos)
bp.add <- c(0, cs)[1:7]
gl <- sum(max.pos)
c.col <- c("darkblue","dodgerblue1","darkblue", "dodgerblue1", "darkblue","dodgerblue1","darkblue", "dodgerblue1")
cp.start <- bp.add[1:7]+(0.5*max.pos)

### var.label is the name of the variable for labeling plot y axis
### win.var is a data frame of chr pos value pvalue
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
        points(c.pos, c.dat,type="p", col=c.col[up.c])
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


pdf("fst.differences.with.pvals.pdf", width=12, height=14)
par(mfcol=c(4,1))
pval.cutoff=-log10(0.05/nrow(fst))
win.var <- cbind(fst[,2:3], bw.diff, -log10(bw.kw))
genome.plot.pval(var.label="between-within fst difference", win.var=win.var, plot.name="between-within", pval.cutoff=pval.cutoff)	

win.var <- cbind(fst[,2:3], fp.diff, -log10(fp.kw))
genome.plot.pval(var.label="formosa-pubescens fst difference", win.var=win.var, plot.name="formosa-pubescens", pval.cutoff=pval.cutoff)	

win.var <- cbind(fst[,2:3], bf.diff, -log10(bf.kw))
genome.plot.pval(var.label="between-formosa fst difference", win.var=win.var, plot.name="between-formosa", pval.cutoff=pval.cutoff)	
	
win.var <- cbind(fst[,2:3], bp.diff, -log10(bp.kw))
genome.plot.pval(var.label="between-pubescens fst difference", win.var=win.var, plot.name="between-pubescens", pval.cutoff=pval.cutoff)	
dev.off()


#### very little difference between bf and bp comparisons.

### make a summary table of all this part and output

sum.fst <- cbind(fst[,1:4], all.fst, all.kw[5:8])
colnames(sum.fst)[5:8]<- paste("fst.",colnames(sum.fst)[5:8], sep="")
colnames(sum.fst)[9:12]<- paste("pval.",colnames(sum.fst)[9:12], sep="")
write.table(sum.fst, file="fst.summary.by.gene.txt", quote=FALSE, row.names=FALSE)



	
#######################################################
### step three: permute
#######################################################



#######################################################
### step four: recombination rate relationships?
#######################################################


