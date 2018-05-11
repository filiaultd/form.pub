### small description of scripts in this folder
### input files, output files
### if two scripts given, one is the submission script to run the other on the cluster



00 -> helper R functions for calculating divergence and pi
01/02 -> from sync file, get polymoprphic sites (generates no169.pools.poly.sync)
04 -> helper R functions for 25 (originally developed for 03)

###############################################################
### Script 03, 05-08 are obsolete but here is a description anyway:
##############################################################

03 -> use filtered sync file to get pool allele freqs  (biallele.freq.Rdata)
05/06 -> from bialleleic allele frequencies, get species means and t-test (biallele.afs.Rdata)
07/08 -> completely obsolete

##############################################################
### From here, started over so could use same data for nucleotide diversity and AFD determination
#############################################################

10/11 -> python script to generate two files from a sync file (in this case, no169.pools.poly.sync)
	no169.pools.poly.sync.goodsites.out -> chr pos 0/1(pass coverage threshold for each pool)
	no169.pools.poly.sync.biallelic.out -> rsync of biallelic SNPs only 

12 -> get individual allele frequencies per pool (file="biallelic.allele.freq.bi.af.Rdata")
	does PCA (../003.output/pca1and2.loadings.allele.freq.pdf, ../003.output/pca.loadings.txt)
	joint SFS (../003.output/joint.SFS.unfolded.pdf, 003.output/SFS.unfolded.form.pub.0.1bin.pdf, etc)
	makes distance matrix (snp.dist.Rdata)

13 -> non existant.  Bad luck??

14/15 -> parses allele frequencies (biallelic.allele.freq.bi.af.Rdata) for significantly-differentiated SNPs.
	kruskal wallace, permutation, t-test (both.mean.af.Rdata)
	does some genome-wide plotting with recomb. rate (003.output/allele.freq.diff.with.recomb.pdf)

16/17 - > calculates divergence (by chromosome and genome-wide) between all pop pairs and then incorporates 18 pi output to calculate Fst (generates 003.output/parsed.all.pairwise.divergence.txt,pairwise.pop.divergence.Rdata)
 
18 -> calculates pi per population (generates ../003.output/all.pop.pi.txt, ../003.output/pop.pi.by.chr.pdf)

19 -> looks at relationship between geo distance and genetic distance.  Also make a map of pool locations
	(generates ../003.output/divergence.by.distance.pdf, 19.geo.analysis.Rdata)

20 -> R functions to do genome rotation in 21 and 22

21 -> overlap between significantly-diverged SNPs and QTL positions

22/23 -> overlap between significantly-diverged SNPs and recombination rate bins
	(rotation.test.snps.by.recomb.bin.pdf)

24 -> shell script (GATK) generates genotypes of Semiaquilegia (from Aquilegia genome paper) at pool bialleleic positions
	(generates SRR.pool.snps.genotype.table.txt, a genotype table)

25 -> gets frequency of derived variants, does joint SFS
	(generates derived.af.Rdata, 003.output/derived.joint.SFS.pdf)

