#!/usr/bin/python
# Merge a VCF to a reference panel
# Switches alleles if strand is unambiguous and:
#	1) chromosomal location matches, but the ref and alt alleles are switched
#	2) chromosomal location matches, but the strand is flipped
# Switches chromosomal position if strand is unambiguous and:
#	1) ID, reference, and alternate allele match
# If alternate allele is missing, and SNP matches by position, and reference allele matches, fills in alternate allele from the reference
# version 1 (09/15/17)


import gzip
import sys
import argparse

class variant:
	def __init__(self):
		chr = '0'
		pos = 0
		ref = 'N'
		alt = 'N'
		id = 'N'

def flipNuc(x):
	y = ""
	for i in x[::-1]:
		if i == "A":
			y = y + "T"
		elif i == "T":
			y = y + "A"
		elif i == "C":
			y = y + "G"
		elif i == "G":
			y = y + "C"
		else:
			y = i
	return(y)

def strandAmbiguous(x,y):
	if (x == "A" and y == "T") or (x == "T" and y == "A") or (x == "G" and y == "C") or (x == "C" and y == "G"):
		return(True)
	else:
		return(False)


def switchGenoVcfLine(x):
	"x is a list of a VCF line with genotypes"
	yl = x
	assert len(yl) >= 9
	gti = yl[8].split(":").index("GT")
	for i in range(9,len(yl)):
		geno = yl[i].split(":")[gti]
		if len(geno) == 3:
			genoA1 = geno[0]
			sep = geno[1]
			genoA2 = geno[2]
			if sep == "|":
				geno = genoA2 + sep + genoA1
			elif geno == "1/1":
				geno = "0/0"
			elif geno == "0/0":
				geno = "1/1"
		elif len(geno) == 1:
			if geno == "0":
				geno = "1"
			elif geno == "1":
				geno = "0"
		else:
			print("Unrecognized genotype in line %s" % s)
			sys.exit(1)
		yl[i] = ":".join(yl[i].split(":")[:gti] + [geno] + yl[i].split(":")[gti+1:])
	return(yl)

parser = argparse.ArgumentParser()
parser.add_argument("--ref", required=True, help="GZipped VCF file of the reference panel to merge")
parser.add_argument("--vcfin", required=True, help="Input VCF to file to match to reference")
parser.add_argument("--vcfout", required=True, help="Output prefix")
args = parser.parse_args()

# ref_vcf_fn = "/u/project/pajukant/malvarez/data/1000G/1000G_vcf/ALL.chr21.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.genotypes.all.vcf.gz"

ref_vcffs = gzip.open(args.ref, 'r')

byID = {}
byPos = {}

print("Storing reference SNPs")
for line in ref_vcffs:
	if line[0] == "#":
		continue
	linel = line.strip().split()
	thisVar = variant()
	thisVar.chr = linel[0]
	thisVar.pos = int(linel[1])
	thisVar.ref = linel[3]
	thisVar.alt = linel[4]
	# Won't handle multi-allelic variants at the moment
	if "," in thisVar.alt:
		continue
	thisVar.id = linel[2]
	chrpos = linel[0] + ":" + linel[1]
	byID[thisVar.id] = thisVar # Variants are unique
	if chrpos in byPos:
		byPos[chrpos].append(thisVar)
	else:
		byPos[chrpos] = [thisVar]

print("Finished storing reference SNPs")
test_vcf_fn = args.vcfin
test_vcf_fs = gzip.open(test_vcf_fn, 'r')
fixd_vcf_fn = args.vcfout + ".fixed.vcf"
fixd_vcf_fs = open(fixd_vcf_fn, 'w')
unfixd_vcf_fn = args.vcfout + ".unfixed.vcf"
unfixd_vcf_fs = open(unfixd_vcf_fn, 'w')

# Will output corrected SNPs, and SNPs that could not be corrected
# This involves changing IDs and positions, and flipping strand where appropriate
# See if SNP is in reference by chromosomal position, or by ID
# If it is any of the above, see if the SNP matches by allele.
# If it matches by allele, change the chromosomal position or ID accordingly.
# Also flip or switch allele if necessary

nFound = 0
nNoMatch = 0
nAlleleMatch = 0
nFlpd = 0
nSwtd = 0
nStrndAmbg = 0
nLocSwtd = 0
for line in test_vcf_fs:
	if line[0] == "#":
		fixd_vcf_fs.write(line)
		unfixd_vcf_fs.write(line)
		continue
	linel = line.strip().split()
	fixd = False
	ref = linel[3]
	alt = linel[4]
	id = linel[2]
	chrpos = linel[0] + ":" + linel[1]
	refalt = ref + alt
	refaltFlp = flipNuc(ref) + flipNuc(alt)
	refaltSwt = alt + ref
	flipped = False
	switched = False
	IDMatch = False
	chrposMatch = False
	strandAmbg = False
	found = False
	if chrpos in byPos:
		chrposMatch = True
	if id in byID:
		IDMatch = True
	if chrpos in byPos:
		mi = -1
		ai = -1
		# If there are multiple reference panel SNPs in this position, loop through them to find a match and get the index
		# Can match by ID or by allele
		for v in range(len(byPos[chrpos])):
			if byPos[chrpos][v].id == id:
				mi = v
			qrefalt = byPos[chrpos][v].ref + byPos[chrpos][v].alt
			# If alternate allele is missing, fill in if the reference allele matches any allele from the panel
			if alt == ".":
				if ref == byPos[chrpos][v].ref:
					alt = byPos[chrpos][v].alt
				elif ref == byPos[chrpos][v].alt:
					alt = byPos[chrpos][v].ref
				linel[4] = alt
				refalt = ref + alt
				refaltFlp = flipNuc(ref) + flipNuc(alt)
				refaltSwt = alt + ref
			strandAmbg = strandAmbiguous(ref, alt)
			if strandAmbg:
				break
			if refalt == qrefalt:
				ai = v
				break
			elif refaltFlp == qrefalt:
				ai = v
				flipped = True
				break
			elif refaltSwt == qrefalt:
				ai = v
				switched = True
				break
		if ai > -1:
			if linel[2] != byPos[chrpos][ai].id:
				IDSwtchd = True
				linel[2] = byPos[chrpos][ai].id
			if flipped and not switched:
				linel[3] = flipNuc(linel[3])
				linel[4] = flipNuc(linel[4])
			elif switched and not flipped:
				a2 = linel[4]
				linel[4] = linel[3]
				linel[3] = a2
				linel = switchGenoVcfLine(linel)
			found = True
	# If chrpos doesn't match but ID does
	elif id in byID:
		qrefalt = byID[id].ref + byID[id].alt
		qchrpos = byID[id].chr + ":" + str(byID[id].pos)
		strandAmbg = strandAmbiguous(ref, alt)
		if not strandAmbg and refalt == qrefalt and chrpos == qchrpos:
			found = True
			nLocSwtd = nLocSwtd + 1
			linel[0] = byID[id].chr
			linel[1] = byID[id].pos
	
	if strandAmbg:
		nStrndAmbg = nStrndAmbg + 1
	elif found and not flipped and not switched:
		nAlleleMatch = nAlleleMatch + 1
	elif found and flipped:
		nFlpd = nFlpd + 1
	elif found and switched:
		nSwtd = nSwtd + 1
	
	if not found:
		nNoMatch = nNoMatch + 1
		unfixd_vcf_fs.write("\t".join(linel) + "\n")
	else:
		nFound = nFound + 1
		fixd_vcf_fs.write("\t".join(linel) + "\n")

unfixd_vcf_fs.close()
fixd_vcf_fs.close()

print("Variants matched with reference:\t%s" % str(nFound))
print("Variants matching both ref and alt alleles:\t%s" % str(nAlleleMatch))
print("Variants matching with ref and alt switched:\t%s" % str(nSwtd))
print("Variants matching with strand flipped:\t%s" % str(nFlpd))
print("Variants not matched with reference:\t%s" % str(nNoMatch))
print("Variants with strand ambiguous alleles:\t%s" % str(nStrndAmbg))
