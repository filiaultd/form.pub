### functions for calculating divergence
### DLF 23April18

############ functions for divergence calcs in
############ 16.pairwise.pop.divergence.R
############ can also be co-opted to do in windows or by chromosome!

### n1, n2 is number of chromosomes per pop (n.indiv*2)
### p1, p2 is proportion of allele in pop

get.mismatch.no <- function(p1, p2, n1,n2){
		ncombos <- n1*n2
		nsamep <- ncombos*p1*p2
		nsameq <- ncombos*(1-p1)*(1-p2)
		ndif <- ncombos-(nsamep+nsameq)
		return(ndif)	
		}

### bi.af is frequency of biallelic snps
### n1, n2 is number of chromosomes per pop (n.indiv*2)
### p1, p2 is proportion of allele in pop
### gsites is the matrix of good sites and their "pass" status
### bi.snps is the data frame of biallelic snps (gives the pos for bi.af)
	
calc.div <- function(bi.af, pop1,pop2, n1, n2, gsites, bi.snps){
	## first get subset of good sites for this comparison
	gsites.pats <- unique(gsites[,3])
	gs.pats.up <- sapply(gsites.pats, function(x){
		pat.g <- substr(x,pop1,pop1)==1 && substr(x,pop2,pop2)==1
		return(pat.g)
	})
	gs.pats.up <- gsites.pats[gs.pats.up]
	gsites.up <- gsites[gsites[,3]%in%gs.pats.up,1:2]
	colnames(gsites.up) <- c("chr","pos")
	
	## get allele freq data subset
	## also need to subset bi.af for only good positions
	## positions are same as bi.snps positions	
	up.af <- bi.af[,c(pop1,pop2)]  ### need to add positions?
	up.af <- cbind(bi.snps[,1:2], up.af)
	good.af <- up.af[up.af[,3]!=up.af[,4],] ###get rid of nonpoly and NA data
	good.af <- good.af[is.na(good.af[,3])==FALSE,]
	colnames(good.af) <- c("chr","pos","af1","af2")
	
	gaf.s <- split(good.af, good.af$chr)
	gsites.s <- split(gsites.up, gsites.up$chr)
	for (up in 1:length(gaf.s)){
		up.gaf <- gaf.s[[up]]
		up.gsites <- gsites.s[[up]]
		g.up.gaf <- up.gaf[up.gaf$pos%in%up.gsites$pos,]
		gaf.s[[up]] <- g.up.gaf
	}
	good.af <- do.call(rbind, gaf.s)
	good.pos <- good.af[,1:2]
	good.af <- good.af[,3:4]
	
	### get total number of mismatch comparisons
	out.mismatch <- rep(NA, nrow(good.af))
	for(up.snp in 1:nrow(good.af)){
		up.mm <- get.mismatch.no(p1=good.af[up.snp,1], p2=good.af[up.snp,2], n1=n1,n2=n2)
	out.mismatch[up.snp]<-up.mm
	}
	total.mm <- sum(out.mismatch)
	out.mismatch <- cbind(good.pos, out.mismatch)
	chr.mm <- aggregate(out.mismatch~chr, data=out.mismatch, sum)
	
	### now need number of good sites and number of total comparisons
	n.good.sites <- unlist(lapply(gsites.s,nrow))
	chr.comps <- n.good.sites*(n1*n2)
	total.comps <- sum(chr.comps)

	### calculate divergence
	total.div <- total.mm/total.comps
	chr.div <- chr.mm[,2]/chr.comps
	all.div <- c(chr.div, total.div)
	return(all.div)
	}
	
##################################################################	
############ functions for pi calcs in
############ 18.nucleotide.diversity.within.R
############ can also be co-opted to do in windows or by chromosome!
#######################################################################


### bi.af is frequency of biallelic snps
### n1 is number of chromosomes per pop (n.indiv*2)
### p1 is proportion of allele in pop
### gsites is the matrix of good sites and their "pass" status
### bi.snps is the data frame of biallelic snps (gives the pos for bi.af)


calc.pi <- function(bi.af, pop1, n1, gsites, bi.snps){

## first get subset of good sites for this comparison
	gsites.pats <- unique(gsites[,3])
	gs.pats.up <- sapply(gsites.pats, function(x){
		pat.g <- substr(x,pop1,pop1)==1
		return(pat.g)
	})
	gs.pats.up <- gsites.pats[gs.pats.up]
	gsites.up <- gsites[gsites[,3]%in%gs.pats.up,1:2]
	colnames(gsites.up) <- c("chr","pos")
	
	## get allele freq data subset
	## also need to subset bi.af for only good positions
	## positions are same as bi.snps positions	
	up.af <- bi.af[,pop1]  ### need to add positions?
	up.af <- cbind(bi.snps[,1:2], up.af)
	good.af <- up.af[up.af[,3]%in%c(0,1)==FALSE,] ###get rid of nonpoly and NA data
	good.af <- good.af[is.na(good.af[,3])==FALSE,]
	colnames(good.af) <- c("chr","pos","af1")
	
	gaf.s <- split(good.af, good.af$chr)
	gsites.s <- split(gsites.up, gsites.up$chr)
	for (up in 1:length(gaf.s)){
		up.gaf <- gaf.s[[up]]
		up.gsites <- gsites.s[[up]]
		g.up.gaf <- up.gaf[up.gaf$pos%in%up.gsites$pos,]
		gaf.s[[up]] <- g.up.gaf
	}
	good.af <- do.call(rbind, gaf.s)
	good.pos <- good.af[,1:2]
	good.af <- good.af[,3]

### get total number of mismatch comparisons
	out.mismatch <- rep(NA, length(good.af))
	for(up.snp in 1:length(good.af)){
		up.mm <- pi.get.mismatch.no(p1=good.af[up.snp], n1=n1)
	out.mismatch[up.snp]<-up.mm
	}
	total.mm <- sum(out.mismatch)
	out.mismatch <- cbind(good.pos, out.mismatch)
	chr.mm <- aggregate(out.mismatch~chr, data=out.mismatch, sum)
	
	### now need number of good sites and number of total comparisons
	n.good.sites <- unlist(lapply(gsites.s,nrow))
	chr.comps <- n.good.sites*(choose(n1,2))
	total.comps <- sum(chr.comps)

	### calculate divergence
	total.pi <- total.mm/total.comps
	chr.pi <- chr.mm[,2]/chr.comps
	all.pi <- c(chr.pi, total.pi)
	return(all.pi)
	}


pi.get.mismatch.no <- function(p1,n1){
		ncombos <- choose(n1, 2)
		ndiff <- (2)*(ncombos*(1-p1)*(p1)) ## since all biallelic
		return(ndiff)	
		}
