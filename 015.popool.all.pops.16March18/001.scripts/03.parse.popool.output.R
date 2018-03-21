### parse output of popool (pi and theta)
### in windows, by gene
### DLF 20MAR18

setwd("/Volumes/aquilegia/004.pooled.sequencing/015.popool.all.pops.16March18/003.output/")
index.file <- ("/Volumes/aquilegia/004.pooled.sequencing/001.data/001.bam.formatted.files/READ.ME.for.Aquilegia.filenames.csv")
gtf.file <- ("/Volumes/aquilegia/001.common.reference.files/018.Acoerulea_v3.1_annotation_20150930/Acoerulea_396_v3.1_genes_only.gtf")


### prepare index so know which pool is which species
pool.index <- read.csv(index.file, stringsAsFactors=FALSE)
pool.index <- pool.index[,c(3,4,8)]
pool.index$species <- gsub("\\?","", pool.index$species)
pool.index$pool.name <- gsub("\\?","", pool.index$pool.name)
pool.index <- unique(pool.index)
pool.index$pool.number <- paste("pool",pool.index$pool.number, sep="")
pool.index <- pool.index[pool.index$pool.number%in%c("pool1", "pool6","pool9")==FALSE,]

##########################################################################
### concatenate data and explore
### start with biggest windows first - 100kb windows with 50kb slide
##########################################################################

files <- dir()
files <- files[-grep("params",files)]
up.files <- files[grep("100kb_50kb",files)]



if(exists("out.pi")){rm(out.pi)}
if(exists("out.nsnps")){rm(out.snps)}
for(up in 1:length(up.files)){
	up.file <- up.files[up]
	up.dat <- read.table(up.file,stringsAsFactors=FALSE)
	up.pool <- unlist(strsplit(up.file,split=".", fixed=TRUE))[1]
	colnames(up.dat) <- c("chrom", "pos",paste(up.pool,"nsnps", sep="."),"prop",paste(up.pool,"pi", sep="."))

	if(up==1){out.pi <-up.dat[,c(1,2,5)]
		}else {out.pi <- merge(out.pi, up.dat[,c(1,2,5)], all.x=TRUE)}
	
	if(up==1){out.nsnps <- up.dat[,c(1,2,3)]
		}else {out.nsnps <- merge(out.nsnps, up.dat[,c(1,2,3)], all.x=TRUE)}
		
	}
	
pi <- out.pi[grep("Chr",out.pi$chr),]
colnames(pi) <- gsub(".pi","",colnames(pi))
for(up in 3:12){class(pi[,up])<-"numeric"}
	
#### look at histograms of pi values per pool
pdf("../004.parsed.data/pi.by.pool.100kb_50kb.pdf", width=16, height=7)
par(mfrow=c(2,5))
for(up in 3:12){
	hist(pi[,up], main=colnames(pi)[up])
	mtext(side=3, pool.index[match(colnames(pi)[up], pool.index$pool.number),3],cex=0.8)
}
dev.off()
### pretty similar, some skewed slightly lower

#### distribution by chromosome
pdf("../004.parsed.data/pi.by.chr.100kb_50kb.pdf", width=16, height=7)
par(mfrow=c(2,5))
for(up in 3:12){
	boxplot(as.numeric(pi[,up])~pi$chr, main=colnames(pi)[up])
}
dev.off()
### definite elevation on chromosome four in all pools

#### means by species, species differences?
p.pool <- pool.index[pool.index$species=="pubescens",2]
p.dat <- pi[,colnames(pi)%in%p.pool]
p.mean <- apply(p.dat,1,mean)
f.pool <- pool.index[pool.index$species=="formosa",2]
f.dat <- pi[,colnames(pi)%in%f.pool]
f.mean <- apply(f.dat,1,mean)
scatterhist(f.mean, p.mean, xlab="mean pi formosa", ylab="mean pi pubescens")

pi.ttest <- apply(pi,1, function(x){
	p.dat <- as.numeric(x[colnames(pi)%in%p.pool])
	f.dat <- as.numeric(x[colnames(pi)%in%f.pool])
	ttt <- try(t.test(p.dat, f.dat))
	if(length(ttt)==1){t.up <- NA
	} else {
	t.up <- t.test(p.dat, f.dat, na.rm=TRUE)$p.value
	}
})

pi$pval <- pi.ttest
pi$p.mean <- p.mean
pi$f.mean <- f.mean
pi$diff <- pi$p.mean-pi$f.mean
hist(pi$diff, xlab="pub pi - form pi")  ### a bit skewed towards pub having higher pi (but very little)
pi.sig <- pi[pi$pval<=(0.05),]
pi.sig <- pi.sig[is.na(pi.sig[,1])==FALSE,]
hist(pi.sig$diff, xlab="pub pi - form pi")  ### again, more skewed towards pub.  848 windows, but not multiple testing corrected...


#### plot pi across genome
#### overplot all?  or mean?

pi.s <- split(pi, pi$chrom)
win.var <- pi.s

max.pos <- lapply(win.var, function(x){max(x[,2])})
max.pos <- unlist(max.pos)
cs <- cumsum(max.pos)
bp.add <- c(0, cs)
gl <- sum(max.pos)
c.col <- c("darkblue","dodgerblue1","darkblue", "dodgerblue1", "darkblue","dodgerblue1","darkblue", "dodgerblue1")
cp.start <- bp.add[1:7]+(0.5*max.pos)
up.pos <- lapply(win.var, function(x)x[,2])

pdf(file="../004.parsed.data/pi.100kb.50kb.across.genome.pdf", width=8.5, height=12)
par(mfcol=c(3,1))
for (up.sp in 14:16){
        par(mar=c(2.5,4.5,1.5,0.5))
        up.sn <- colnames(pi)[up.sp]
        up.prefix <- colnames(pi)[up.sp]
        up.dat <- lapply(win.var,function(x){x[,colnames(x)==up.prefix]})
        up.dat <- lapply(up.dat, function(x){x*100})
        max.pi <- max(unlist(up.dat), na.rm=TRUE)
        min.pi <- min(unlist(up.dat), na.rm=TRUE)

        #max.pi <- 1
        plot(0,0, xlim=c(0,gl),ylim=c(min.pi,max.pi),type="n", xaxt="n",ylab="",xlab="")
        for(up.c in 1:7){
                c.dat <- up.dat[[up.c]]
                c.bp.add <- bp.add[up.c]
                c.pos <- up.pos[[up.c]] + c.bp.add
                points(c.pos, c.dat,type="p", col=c.col[up.c])
        }
        abline(v=c(0,cs), col="grey50",lty=3)
        mtext(up.sn, side=3, line=-1.5, cex=0.8, adj=0.975, font=4)
        mtext(1:7, 1,at=cp.start, line=0.25, cex=0.7)
        mtext("chromosome", 1, line=1.3, cex=0.7)
        mtext("percentage pairwise differences",2,las=3,cex=0.7, line=2.2)
}
dev.off()

#### save this variable
 kb100.pi <- pi
 save(kb100.pi, file="../004.parsed.data/kb100.pi.Rdata")




##########################################################################
### concatenate data and explore
### repeat for smaller windows - 10kb windows with 5kb slide
##########################################################################

up.files <- files[grep("10kb_5kb",files)]



if(exists("out.pi")){rm(out.pi)}
if(exists("out.nsnps")){rm(out.nsnps)}
for(up in 1:length(up.files)){
	up.file <- up.files[up]
	up.dat <- read.table(up.file,stringsAsFactors=FALSE)
	up.pool <- unlist(strsplit(up.file,split=".", fixed=TRUE))[1]
	colnames(up.dat) <- c("chrom", "pos",paste(up.pool,"nsnps", sep="."),"prop",paste(up.pool,"pi", sep="."))

	if(up==1){out.pi <-up.dat[,c(1,2,5)]
		}else {out.pi <- merge(out.pi, up.dat[,c(1,2,5)], all.x=TRUE)}
	
	if(up==1){out.nsnps <- up.dat[,c(1,2,3)]
		}else {out.nsnps <- merge(out.nsnps, up.dat[,c(1,2,3)], all.x=TRUE)}
		
	}
	
pi <- out.pi[grep("Chr",out.pi$chr),]
colnames(pi) <- gsub(".pi","",colnames(pi))
for(up in 3:12){class(pi[,up])<-"numeric"}

	
#### look at histograms of pi values per pool
pdf("../004.parsed.data/pi.by.pool.10kb_5kb.pdf", width=16, height=7)
par(mfrow=c(2,5))
for(up in 3:12){
	hist(pi[,up], main=colnames(pi)[up])
	mtext(side=3, pool.index[match(colnames(pi)[up], pool.index$pool.number),3],cex=0.8)
}
dev.off()
### pretty similar, some skewed slightly lower

#### distribution by chromosome
pdf("../004.parsed.data/pi.by.chr.10kb_5kb.pdf", width=16, height=7)
par(mfrow=c(2,5))
for(up in 3:12){
	boxplot(as.numeric(pi[,up])~pi$chr, main=colnames(pi)[up])
}
dev.off()
### definite elevation on chromosome four in all pools


#### means by species, species differences?
p.pool <- pool.index[pool.index$species=="pubescens",2]
p.dat <- pi[,colnames(pi)%in%p.pool]
p.mean <- apply(p.dat,1,mean)
f.pool <- pool.index[pool.index$species=="formosa",2]
f.dat <- pi[,colnames(pi)%in%f.pool]
f.mean <- apply(f.dat,1,mean)
scatterhist(f.mean, p.mean, xlab="mean pi formosa", ylab="mean pi pubescens") ##  shift higher in pub?

pi.ttest <- apply(pi,1, function(x){
	p.dat <- as.numeric(x[colnames(pi)%in%p.pool])
	f.dat <- as.numeric(x[colnames(pi)%in%f.pool])
	ttt <- try(t.test(p.dat, f.dat))
	if(length(ttt)==1){t.up <- NA
	} else {
	t.up <- t.test(p.dat, f.dat, na.rm=TRUE)$p.value
	}
})

pi$pval <- pi.ttest
pi$p.mean <- p.mean
pi$f.mean <- f.mean
pi$diff <- pi$p.mean-pi$f.mean
hist(pi$diff, xlab="pub pi - form pi")  ### a bit skewed towards pub having higher pi (but very little)
pi.sig <- pi[pi$pval<=(0.05),]
pi.sig <- pi.sig[is.na(pi.sig[,1])==FALSE,]
hist(pi.sig$diff, xlab="pub pi - form pi")  ### again, more skewed towards pub.  848 windows, but not multiple testing corrected...

#### plot pi across genome
#### plotting mean and differences

pi.s <- split(pi, pi$chrom)
win.var <- pi.s

max.pos <- lapply(win.var, function(x){max(x[,2])})
max.pos <- unlist(max.pos)
cs <- cumsum(max.pos)
bp.add <- c(0, cs)
gl <- sum(max.pos)
c.col <- c("darkblue","dodgerblue1","darkblue", "dodgerblue1", "darkblue","dodgerblue1","darkblue", "dodgerblue1")
cp.start <- bp.add[1:7]+(0.5*max.pos)
up.pos <- lapply(win.var, function(x)x[,2])

pdf(file="../004.parsed.data/pi.10kb.5kb.across.genome.pdf", width=8.5, height=12)
par(mfcol=c(3,1))
for (up.sp in 14:16){
        par(mar=c(2.5,4.5,1.5,0.5))
        up.sn <- colnames(pi)[up.sp]
        up.prefix <- colnames(pi)[up.sp]
        up.dat <- lapply(win.var,function(x){x[,colnames(x)==up.prefix]})
        up.dat <- lapply(up.dat, function(x){x*100})
        max.pi <- max(unlist(up.dat), na.rm=TRUE)
        min.pi <- min(unlist(up.dat), na.rm=TRUE)

        #max.pi <- 1
        plot(0,0, xlim=c(0,gl),ylim=c(min.pi,max.pi),type="n", xaxt="n",ylab="",xlab="")
        for(up.c in 1:7){
                c.dat <- up.dat[[up.c]]
                c.bp.add <- bp.add[up.c]
                c.pos <- up.pos[[up.c]] + c.bp.add
                points(c.pos, c.dat,type="p", col=c.col[up.c])
        }
        abline(v=c(0,cs), col="grey50",lty=3)
        mtext(up.sn, side=3, line=-1.5, cex=0.8, adj=0.975, font=4)
        mtext(1:7, 1,at=cp.start, line=0.25, cex=0.7)
        mtext("chromosome", 1, line=1.3, cex=0.7)
        mtext("percentage pairwise differences",2,las=3,cex=0.7, line=2.2)
}
dev.off()
### crazy noisy!!!


#### save this variable
 kb10.pi <- pi
 save(kb10.pi, file="../004.parsed.data/kb10.pi.Rdata")


##########################################################################
### concatenate data and explore
### Finally look at values by gene
##########################################################################

up.files <- files[grep("gene",files)]
up.files <- up.files[grep("pi",up.files)]



if(exists("out.pi")){rm(out.pi)}
if(exists("out.nsnps")){rm(out.nsnps)}
for(up in 1:length(up.files)){
	up.file <- up.files[up]
	print(up.file)
	up.dat <- read.table(up.file,stringsAsFactors=FALSE)
	up.pool <- unlist(strsplit(up.file,split=".", fixed=TRUE))[2]
	colnames(up.dat) <- c("gene",paste(up.pool,"nsnps", sep="."),"prop",paste(up.pool,"pi", sep="."))

	if(up==1){out.pi <-up.dat[,c(1,4)]
		}else {out.pi <- merge(out.pi, up.dat[,c(1,4)], all.x=TRUE)}
	
	if(up==1){out.nsnps <- up.dat[,c(1,2)]
		}else {out.nsnps <- merge(out.nsnps, up.dat[,c(1,2)], all.x=TRUE)}
		
	}
	
pi <- out.pi[grep("G",out.pi$gene),]
colnames(pi) <- gsub(".pi","",colnames(pi))
for(up in 2:11){class(pi[,up])<-"numeric"}

### parse to get chromosome and positions

gtf <- read.table(gtf.file, stringsAsFactors=FALSE)
gtf <- gtf[,c(1,4,5,10)]
colnames(gtf) <- c("chrom", "start", "stop", "gene")
pi <- merge(gtf,pi, by="gene",all=TRUE)
pi <- pi[grep("Chr", pi$chrom),]  #19550 primary gene models


#### look at histograms of pi values per pool
pdf("../004.parsed.data/pi.by.pool.genes.pdf", width=16, height=7)
par(mfrow=c(2,5))
for(up in 5:14){
	hist(pi[,up], main=colnames(pi)[up])
	mtext(side=3, pool.index[match(colnames(pi)[up], pool.index$pool.number),3],cex=0.8)
}
dev.off()
### pretty similar, some skewed slightly lower

#### distribution by chromosome
pdf("../004.parsed.data/pi.by.chr.genes.pdf", width=16, height=7)
par(mfrow=c(2,5))
for(up in 5:14){
	boxplot(as.numeric(pi[,up])~pi$chr, main=colnames(pi)[up])
	mtext(side=3, pool.index[match(colnames(pi)[up], pool.index$pool.number),3],cex=0.8)
}
dev.off()
### definite elevation on chromosome four in all pools


#### means by species, species differences?
p.pool <- pool.index[pool.index$species=="pubescens",2]
p.dat <- pi[,colnames(pi)%in%p.pool]
p.mean <- apply(p.dat,1,mean)
f.pool <- pool.index[pool.index$species=="formosa",2]
f.dat <- pi[,colnames(pi)%in%f.pool]
f.mean <- apply(f.dat,1,mean)
scatterhist(f.mean, p.mean, xlab="mean pi formosa", ylab="mean pi pubescens") 

pi.ttest <- apply(pi,1, function(x){
	p.dat <- as.numeric(x[colnames(pi)%in%p.pool])
	f.dat <- as.numeric(x[colnames(pi)%in%f.pool])
	ttt <- try(t.test(p.dat, f.dat))
	if(length(ttt)==1){t.up <- NA
	} else {
	t.up <- t.test(p.dat, f.dat, na.rm=TRUE)$p.value
	}
})

pi$pval <- pi.ttest
pi$p.mean <- p.mean
pi$f.mean <- f.mean
pi$diff <- pi$p.mean-pi$f.mean
hist(pi$diff, xlab="pub pi - form pi")  ### a bit skewed towards pub having higher pi (but very little)
pi.sig <- pi[pi$pval<=(0.05),]
pi.sig <- pi.sig[is.na(pi.sig[,1])==FALSE,]
hist(pi.sig$diff, xlab="pub pi - form pi")  ### again, more skewed towards pub.  5712 genes, but not multiple testing corrected...

pdf(file="../004.parsed.data/pval.vs.diff.genes.pdf", width=8, height=8)
plot(-log10(pi.sig$pval), pi.sig$diff, xlab="-log10 Pval", ylab="mean pi difference (pub-form)")
dev.off()
### so, when there are significant differences in pi between species, these tend to be cases in which pi is higher in pubescens.




#### plot pi across genome
#### plotting mean and differences

pi.s <- split(pi, pi$chrom)
win.var <- pi.s

max.pos <- lapply(win.var, function(x){max(x[,4])})
max.pos <- unlist(max.pos)
cs <- cumsum(max.pos)
bp.add <- c(0, cs)
gl <- sum(max.pos)
c.col <- c("darkblue","dodgerblue1","darkblue", "dodgerblue1", "darkblue","dodgerblue1","darkblue", "dodgerblue1")
cp.start <- bp.add[1:7]+(0.5*max.pos)
up.pos <- lapply(win.var, function(x)x[,3])

pdf(file="../004.parsed.data/pi.genes.across.genome.pdf", width=8.5, height=12)
par(mfcol=c(3,1))
for (up.sp in 16:18){
        par(mar=c(2.5,4.5,1.5,0.5))
        up.sn <- colnames(pi)[up.sp]
        up.prefix <- colnames(pi)[up.sp]
        up.dat <- lapply(win.var,function(x){x[,colnames(x)==up.prefix]})
        up.dat <- lapply(up.dat, function(x){x*100})
        max.pi <- max(unlist(up.dat), na.rm=TRUE)
        min.pi <- min(unlist(up.dat), na.rm=TRUE)

        #max.pi <- 1
        plot(0,0, xlim=c(0,gl),ylim=c(min.pi,max.pi),type="n", xaxt="n",ylab="",xlab="")
        for(up.c in 1:7){
                c.dat <- up.dat[[up.c]]
                c.bp.add <- bp.add[up.c]
                c.pos <- up.pos[[up.c]] + c.bp.add
                points(c.pos, c.dat,type="p", col=c.col[up.c])
        }
        abline(v=c(0,cs), col="grey50",lty=3)
        mtext(up.sn, side=3, line=-1.5, cex=0.8, adj=0.975, font=4)
        mtext(1:7, 1,at=cp.start, line=0.25, cex=0.7)
        mtext("chromosome", 1, line=1.3, cex=0.7)
        mtext("percentage pairwise differences",2,las=3,cex=0.7, line=2.2)
}
dev.off()
### crazy noisy!!!  but no obvious patterns now like seen in the windows...
### pools 12 and 13 aren't complete - need to rerun...































################################# helper functions - run before script ##########################
# from https://www.r-bloggers.com/example-8-41-scatterplot-with-marginal-histograms/
# slight mod to fix placement of x and y labels
scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(4,4,1,1))
  plot(x,y,xlab=xlab, ylab=ylab)
  par(mar=c(0,4,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(4,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(4,4,0,0))
  #mtext(xlab, side=1, line=2, outer=TRUE, adj=0.4)
  #mtext(ylab, side=2, line=2, outer=TRUE, adj=0.4)

}
