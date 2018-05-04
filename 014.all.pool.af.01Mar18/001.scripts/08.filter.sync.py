import sys
import os

def main():
	filepath=sys.argv[1] ### whole genome sync file path
	#outfile="test.out"
	outfile=sys.argv[2]  ### output file name/path
	#minCov=int(sys.argv[3])  ### minimum coverage for a pool to be considered
	#minCov=4

	if not os.path.isfile(filepath):
		print("File path {} does not exist.  Exiting".format(filepath))
		sys.exit()
	
	fileOpen = open(outfile, "w")

	poly_sites={}
	with open(filepath) as upfile:
		for lineOut in upfile:
			line=lineOut.rstrip('\n')
			#print("line {} is up next".format(line))
			line_s=line.split("\t")
			#print("splitline is {}".format(line_s))
			ref=line_s[2]
			#print("ref allele is {}".format(ref))			

			### cov of each pool
			alleles = line_s[3:13]
			#print("alleles are {}".format(alleles))
			pCov=[]
			for up_allele in alleles:
				poolCov(up_allele, pCov)			

			### generate output line
			pos=line_s[0:1]
			pCovOut = pos + pCov
			print(pCovOut)

			fileOpen.write(pCovOut)
	fileOpen.close()

def poolCov(up_allele,pCov):
	ua_s=up_allele.split(":")
        ua_s=map(int, ua_s)
	ua_sum=sum(ua_s)
	pCov.append(ua_sum)
	

if __name__ == '__main__':
	main()
