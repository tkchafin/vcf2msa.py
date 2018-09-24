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
		name = contig[0].split()[0]
		reference[name] = contig[1] #don't forget: 0-based index here, 1-based in VCF

	vfh = vcf.Reader(open(params.vcf, 'r'))

	#grab contig sizes
	contigs = dict()
	for c,s in vfh.contigs.items():
		contigs[s.id] = s.length

	#Data structure with:
	#For each chromosome:
	#	for each position:
	#		mpileup depth
	#		VCF ALT allele
	#Dict of lists [N] of lists [2]

	data = dict()
	for chrom, seq in reference.items():
		l1 = [None] * len(seq)
		l2 = [None] * len(seq)
		data[chrom] = [l1, l2]

	#First parse the mpileup file
	with open(params.pileup, 'r') as PILEUP:
		try:
			for l in PILEUP:
				line = l.split()
				chrom = line[0]
				pos = int(line[1])-1
				depth = line[3]
				if chrom not in data:
					pass
				else:
					if pos < len(data[chrom][0]):
						data[chrom][0][pos] = depth
		except IOError as e:
			print("Could not read file %s: %s"%(params.pileup,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(params.pileup,e))
			sys.exit(1)
		finally:
			PILEUP.close()

	#Parse VCF to grab genotypes
	#Note this doesn't pay attention to samples. just keeps all ALT alleles
	for rec in vfh:
		chrom = rec.CHROM
		pos = int(rec.POS)-1
		allele = rec.ALT
		chosen = rec.ALT[0]

		#if more than one ALT allele, choose shortest
		#(i.e. SNP over indel, or shortest indel)
		#or keep first if equal
		if len(allele) > 1:
			for a in allele:
				if len(a) < len(chosen):
					chosen = a

		if chrom not in data:
			pass
		else:
			if pos < len(data[chrom][1]):
				data[chrom][1][pos] = chosen

	#For each position in chromosome, write new one for sample
	for chrom, seq in reference.items():
		print(">%s_%s"%(chrom, params.sample))
		for position, nuc in enumerate(seq):
			if int(data[chrom][0][position]) >= params.cov:
				if data[chrom][1][position] != None:
					if data[chrom][1][position] == "*":
						pass
					else:
						print(str(data[chrom][1][position]), end="")
				else:
					print(nuc, end="")
			else:
				print("N",end="")
		print()


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

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'r:v:m:c:R:hs:', \
			["vcf=", "help", "ref=", "mpileup=","cov=","region=","sample="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.ref=None
		self.pileup=None
		self.cov=0
		self.sample="sample"

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
			elif opt in ('m','mpileup'):
				self.pileup=arg
			elif opt in ('r','ref'):
				self.ref=arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('s','sample'):
				self.sample=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Must provide VCF file <-v,--vcf>")
		if not self.ref:
			self.display_help("Must provide reference FASTA file <-r,--ref")
		if not self.pileup:
			self.display_help("Must provide mpilup output <-m,--mpileup>")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\naltRefMaker.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("Description: Builds a reference consensus sequence from VCF and MPILEUP")

		print("""
	Arguments:
		-r,--ref	: Reference genome file (FASTA)
		-v,--vcf	: VCF file containing genotypes
		-m,--mpileup	: Samtools mpileup file for ALL sites (-aa option in mpileup)
		-c,--cov	: Minimum coverage (taken from mpileup) to call, otherwise write "N"
		-s,--sample	: Output sample name (default="sample")
		-h,--help	: Displays help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
