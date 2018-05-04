### calculating divergence between two pools
### prep work is done in scripst 10 and 11 in this directory
### DLF 13April18
### Update 24 April 18 - this script mostly just does site frequency spectra.


input.dir <- "/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input/"
output.dir <- "/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/003.output/"

gs.file <- "no169.pools.poly.sync.goodsites.out"
bi.file <- "no169.pools.poly.sync.biallelic.out"


setwd(input.dir)

gsites <- read.table(gs.file, stringsAsFactors=FALSE, colClasses=c("character","numeric","character"))
## these are "good sites", but still need to be filtered for coverage of individual pools
bi.snps <- read.table(bi.file, stringsAsFactors=FALSE)
## bialleleic snps.  Might also still need to be filtered as gsites
## just check that overlaps gsites after subsetting for 2 pops of interest
pool.index <- read.table("pool.order.index.txt", header=TRUE,stringsAsFactors=FALSE) 
## file that has the number of individuals per pool


pass.pats <- unique(gsites[,3])
pass.pop <- do.call(rbind,sapply(pass.pats,function(x){strsplit(x, split="", fixed=TRUE)}))
### this is the parsed index of all combos of whether a site passed the filters or not.  Will be needed to filter good sites variable

chrs <- unique(bi.snps[,1])[1:7]

## limit to chromosomes only
gsites <- gsites[gsites[,1]%in%chrs,]
bi.snps <- bi.snps[bi.snps[,1]%in%chrs,]

### get allele frequencies
### these are all biallelic, so this is a simpler question

get.af <- function(bi.snps){  ### a data frame of snps
	us.af <- apply(bi.snps, 1, function(s.up){
      	ac <- sapply(s.up[4:length(s.up)], function(dp){out.dp <- strsplit(dp, split=":")})
      	ac <- do.call(rbind, ac)
      	ac <- apply(ac,1, as.numeric)
      	sums <- apply(ac,2, sum)
      	af <- apply(ac,2, function(cu){cu/sum(cu)})
      	af[is.na(af)] <- NA
	  	## check if poly
      	poly <- apply(af,1, function(x){sum(x,na.rm=TRUE)})
      	out.af <- af[poly>0,]
      	if(class(out.af)=="numeric"){out.row <- out.af
      	} else {out.row <- c(out.af[1,])}
      	return(out.row)
      	})
	return(us.af)
}



bi.af <- get.af(bi.snps)
bi.af <- t(bi.af)
save.image("test.Rdata")
#load("test.Rdata")

### save allele frequencies alone so can work with them elsewhere!
### also add positions so can do this by chromosome - that's from bi.snps file

save(bi.af, file="biallelic.allele.freq.bi.af.Rdata")

### PCA of allele frequencies
### won't take NAs/missing data so remove
af.miss <- apply(bi.af, 1, max)
af.complete <- bi.af[is.na(af.miss)==FALSE,]  ##22million positions complete
af.complete <- t(af.complete)
af.pca <- prcomp(af.complete, center = TRUE, scale. = TRUE) 
pca.sum <- summary(af.pca)

pdf("../003.output/pca1and2.loadings.allele.freq.pdf")
plot(af.pca$x[,1:2], col=c("blue","blue","red","blue","blue","red","red","red","red","blue"),pch=19)
text(af.pca$x[,1:2], labels=pool.index$short.name, pos=3)
dev.off()

pdf("../003.output/pca.others.loadings.allele.freq.pdf")
plot(af.pca$x[,3:4], col=c("blue","blue","red","blue","blue","red","red","red","red","blue"),pch=19)
text(af.pca$x[,3:4], labels=pool.index$short.name, pos=3)

plot(af.pca$x[,5:6], col=c("blue","blue","red","blue","blue","red","red","red","red","blue"),pch=19)
text(af.pca$x[,5:6], labels=pool.index$short.name, pos=3)

plot(af.pca$x[,7:8], col=c("blue","blue","red","blue","blue","red","red","red","red","blue"),pch=19)
text(af.pca$x[,7:8], labels=pool.index$short.name, pos=3)

plot(af.pca$x[,9:10], col=c("blue","blue","red","blue","blue","red","red","red","red","blue"),pch=19)
text(af.pca$x[,9:10], labels=pool.index$short.name, pos=3)
dev.off()

pca.loadings <- af.pca$x
write.table(pca.loadings, "../003.output/pca.loadings.txt")

####################################
#### joint site frequency spectrum
####################################

#### mean form/pub

head(bi.af)
head(bi.snps)
which(pool.index$species=="pubescens")
pub.mean.af <- apply(bi.af[,which(pool.index$species=="pubescens")], 1, function(x){mean(x, na.rm=TRUE)})
form.mean.af <- apply(bi.af[,which(pool.index$species=="formosa")], 1, function(x){mean(x, na.rm=TRUE)})

both.mean.af <- cbind(pub.mean.af, form.mean.af)
fold.af <- function(afs){
	for(up in 1:nrow(afs)){
		up.dat <- afs[up,]
		up.max <- max(up.dat, na.rm=TRUE)
		if(up.max>0.5){
			up.dat <- 1-up.dat
			afs[up,]<-up.dat
		}
	}
	return(afs)
}

folded.mean.af <- fold.af(both.mean.af)
colnames(folded.mean.af) <- c("pub.mean","form.mean")

## plot these together

pdf("../003.output/joint.SFS.folded.pdf", width=5, height=5)
smoothScatter(folded.mean.af[,1], folded.mean.af[,2], xlab="formosa mean allele frequency", ylab="pubescens mean allele frequency")
dev.off()


#### each set of pairs

pdf("../003.output/joint.SFS.between.pops.pdf", width=12, height=12)
pairs(bi.af, panel = function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
dev.off()


### get absolute numbers in each bin (or improve the plotting somehow)
### asymmetry between species? need polarization

### also plot unfolded
pdf("../003.output/joint.SFS.unfolded.pdf", width=5, height=5)
smoothScatter(both.mean.af[,1], both.mean.af[,2], xlab="formosa mean allele frequency", ylab="pubescens mean allele frequency")
dev.off()



#####################################################
### simplify this by putting in allele frequency bins
#####################################################


u <- 0.05 #bin size
bin.af <- apply(folded.mean.af, 2, function(x){floor(x/u)*u}) # rounds down to nearest 0.1
colnames(bin.af) <- c("pub","form")
bin.af <- as.data.frame(bin.af)
bin.sfs <- with(bin.af, table(pub,form))

write.table(bin.sfs, file="../003.output/SFS.folded.form.pub.0.05bin.txt")
pdf(file="../003.output/SFS.folded.form.pub.0.05bin.pdf", width=6, height=6)
image(log(bin.sfs), xlab="pub allele freq (log)", ylab="form allele freq(log)")
dev.off()

u <- 0.1 #bin size
bin.af <- apply(folded.mean.af, 2, function(x){floor(x/u)*u}) # rounds down to nearest 0.1
colnames(bin.af) <- c("pub","form")
bin.af <- as.data.frame(bin.af)
bin.sfs <- with(bin.af, table(pub,form))

write.table(bin.sfs, file="../003.output/SFS.folded.form.pub.0.1bin.txt")
pdf(file="../003.output/SFS.folded.form.pub.0.1bin.pdf", width=6, height=6)
image(log(bin.sfs), xlab="pub allele freq (log)", ylab="form allele freq(log)")
dev.off()

#### if one really wants to look at directionality here, one would need to polarize.  Which would take quite a while.
#### can do it if required!

#### get a small accounting of how many SNPs we are talking about at AFD cutoffs
afd.cut <- c(1, 0.95,0.9,0.85)
afd <- abs(folded.mean.af[,1]-folded.mean.af[,2])
afd.cut.no <- sapply(afd.cut, function(x){
	ao <- afd[afd>=x]
	ao <- ao[is.na(ao)==FALSE]
	aol <- length(ao)
	return(aol)
	
})
afd.cut <- cbind(afd.cut, afd.cut.no)
colnames(afd.cut) <- c("cutoff.afd","nsnps")

save.image("12.calc.divergence.RData")

##########################################################################
### simplify this by putting in allele frequency bins - repeat for non-folded
############################################################################
u <- 0.05 #bin size
bin.af <- apply(both.mean.af, 2, function(x){floor(x/u)*u}) # rounds down to nearest 0.1
colnames(bin.af) <- c("pub","form")
bin.af <- as.data.frame(bin.af)
bin.sfs <- with(bin.af, table(pub,form))

write.table(bin.sfs, file="../003.output/SFS.unfolded.form.pub.0.05bin.txt")
pdf(file="../003.output/SFS.unfolded.form.pub.0.05bin.pdf", width=6, height=6)
image(log(bin.sfs), xlab="pub allele freq (log)", ylab="form allele freq(log)")
dev.off()

u <- 0.1 #bin size
bin.af <- apply(both.mean.af, 2, function(x){floor(x/u)*u}) # rounds down to nearest 0.1
colnames(bin.af) <- c("pub","form")
bin.af <- as.data.frame(bin.af)
bin.sfs <- with(bin.af, table(pub,form))

write.table(bin.sfs, file="../003.output/SFS.unfolded.form.pub.0.1bin.txt")
pdf(file="../003.output/SFS.unfolded.form.pub.0.1bin.pdf", width=6, height=6)
image(log(bin.sfs), xlab="pub allele freq (log)", ylab="form allele freq(log)")
dev.off()

#### if one really wants to look at directionality here, one would need to polarize.  Which would take quite a while.
#### can do it if required!




#################################################################
### calculate divergence
### this is now a stand-alone in 16.pairwise.pop.divergence.R
###################################################################

####################################################################
### calculate pi
### also now standalone in 18.pop.pi.R
##################################################################


#####################################################################
### make a distance matrix from allele frequency data
####################################################################

snp.dist <- dist(t(bi.af), upper=TRUE, diag=TRUE)
save(snp.dist, file="snp.dist.Rdata")





