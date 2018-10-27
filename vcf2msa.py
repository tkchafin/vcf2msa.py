#!/usr/bin/python

import re
import sys
import os
import getopt
import vcf

def main():
	params = parseArgs()

	#Grab reference sequence first
	reference = dict()
	for contig in read_fasta(params.ref):
		if params.region:
			if contig[0].split()[0] == params.region.chr:
				name = contig[0].split()[0]
				reference[name] = contig[1] #don't forget: 0-based index here, 1-based in VCF
		else:
			name = contig[0].split()[0]
			reference[name] = contig[1] #don't forget: 0-based index here, 1-based in VCF

	#Get mask sites for each sample
	sampleMask = dict(); #dict of dicts of sets!
	for maskFile in params.mask:
		print(maskFile)
		base=os.path.basename(maskFile)
		samp=base.split(".")[0]
		print(samp)
		with open(maskFile, 'r') as PILEUP:
			try:
				for l in PILEUP:
					line = l.split()
					chrom = line[0]
					pos = int(line[1])-1
					depth = int(line[3])
					#skip if this isn't the targeted chromosome
					if params.region:
						if params.region.chr != chrom:
							continue
						if not (params.region.start <= pos+1 <= params.region.end):
							continue
					if depth < params.cov:
						#Add to sampleMask
						if samp not in sampleMask:
							sampleMask[samp] = dict()
						if chrom not in sampleMask[samp]:
							sampleMask[samp][chrom] = set()
						sampleMask[samp][chrom].add(pos)

			except IOError as e:
				print("Could not read file %s: %s"%(maskFile,e))
				sys.exit(1)
			except Exception as e:
				print("Unexpected error reading file %s: %s"%(maskFile,e))
				sys.exit(1)
			finally:
				PILEUP.close()

	vfh = vcf.Reader(filename=params.vcf)
	samples = list()
	for samp in vfh.samples:
		s = samp.split(".")[0]
		if s not in samples:
			samples.append(s)

	for contig, sequence in reference.items():
		outputs = dict()
		for samp in samples:
			outputs[samp] = str()
		spos = 0
		epos = len(sequence)
		if params.region:
			if contig != params.region.chr:
				continue
			spos = params.region.start-1
			epos = params.region.end-1
			if epos > len(sequence) or spos > (len(sequence)):
				print("WARNING: Specified region outside of bounds: Using full contig",contig)
		#print(spos, " ", epos)

		#loop through each nuc position in contig
		for nuc in range(spos, epos):
			this_pos = dict()
			for samp in samples:
				this_pos[samp] = dict()
			#For each sequence:
			#1. check if any samples have VCF data (write VARIANT)
			for rec in vfh.fetch(contig, nuc+1, nuc+1):
				print(rec)
			#2. check if any samples are masked (write N)
			#3. if NOT masked, write REF



	# #Parse VCF to grab genotypes
	# #Note this doesn't pay attention to samples. just keeps all ALT alleles
	# for rec in vfh:
	# 	chrom = rec.CHROM
	# 	pos = int(rec.POS)-1
	# 	allele = rec.ALT
	# 	chosen = rec.ALT[0]
	#
	# 	#if more than one ALT allele, choose shortest
	# 	#(i.e. SNP over indel, or shortest indel)
	# 	#or keep first if equal
	# 	if len(allele) > 1:
	# 		for a in allele:
	# 			if len(a) < len(chosen):
	# 				chosen = a
	#
	# 	if chrom not in data:
	# 		pass
	# 	else:
	# 		if pos < len(data[chrom][1]):
	# 			data[chrom][1][pos] = chosen
	#
	# #For each position in chromosome, write new one for sample
	# for chrom, seq in reference.items():
	# 	print(">%s_%s"%(chrom, params.sample))
	# 	for position, nuc in enumerate(seq):
	# 		if int(data[chrom][0][position]) >= params.cov:
	# 			if data[chrom][1][position] != None:
	# 				if data[chrom][1][position] == "*":
	# 					pass
	# 				else:
	# 					print(str(data[chrom][1][position]), end="")
	# 			else:
	# 				print(nuc, end="")
	# 		else:
	# 			print("N",end="")
	# 	print()


#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
def read_fasta(fas):

	fh = open(fas)
	try:
		with fh as file_object:
			contig = ""
			seq = ""
			for line in file_object:
				line = line.strip()
				if not line:
					continue
				#print(line)
				if line[0] == ">": #Found a header line
					#If we already loaded a contig, yield that contig and
					#start loading a new one
					if contig:
						yield([contig,seq]) #yield
						contig = "" #reset contig and seq
						seq = ""
					contig = (line.replace(">",""))
				else:
					seq += line
		#Iyield last sequence, if it has both a header and sequence
		if contig and seq:
			yield([contig,seq])
	finally:
		fh.close()

class ChromRegion():
	def __init__(self, s):
		stuff = re.split(':|-', s)
		self.chr = str(stuff[0])
		self.start =  int(stuff[1])
		self.end = int(stuff[2])


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'r:v:m:c:R:hs:', \
			["vcf=", "help", "ref=", "mask=","cov=","reg="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.ref=None
		self.mask = list()
		self.cov=1
		self.region=None

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt in ('v', 'vcf'):
				self.vcf = arg
			elif opt in ('c','cov'):
				self.cov=int(arg)
			elif opt in ('m','mask'):
				self.mask.append(arg)
			elif opt in ('r','ref'):
				self.ref=arg
			elif opt in ('R', 'reg'):
				self.region = ChromRegion(arg)
			elif opt in ('h', 'help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Must provide VCF file <-v,--vcf>")
		if not self.ref:
			self.display_help("Must provide reference FASTA file <-r,--ref")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nvcf2msa.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("Description: Builds multiple sequence alignments from multi-sample VCF")

		print("""
	Arguments:
		-r,--ref	: Reference genome file (FASTA)
		-v,--vcf	: VCF file containing genotypes
			Note: Sample names will be taken as everything before "."
			This is because in my files I have two columns per sample: sample.SNP sample.INDEL
			Feel free to change this if you don't want this behavior
		-m,--mask	: Per-sample mpileup file (one for each sample) for ALL sites (-aa in samtools)
		-c,--cov	: Minimum coverage to call a base as REF [default=1]
		-R,--reg	: Region to sample (e.g. chr1:1-1000)
		-h,--help	: Displays help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
