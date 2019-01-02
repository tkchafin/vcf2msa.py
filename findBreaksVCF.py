#!/usr/bin/python

import re
import sys
import os
import getopt
import vcf

def main():
	params = parseArgs()

	vfh = vcf.Reader(open(params.vcf, 'r'))

	#grab contig sizes
	contigs = dict()
	for c,s in vfh.contigs.items():
		contigs[s.id] = s.length

	regions = list()

	this_chrom = None
	start = int()
	stop = int()
	count = 0
	for rec in vfh:
		if not this_chrom:
			this_chrom = rec.CHROM
			start = 1
			stop = 1
			count = 0
		#If we entered new chromosome, submit old break
		elif this_chrom != rec.CHROM:
			t = tuple([this_chrom, start, contigs[this_chrom]])
			regions.append(t)
			this_chrom = rec.CHROM
			start = 1
			stop = 1
			count = 0

		#if this SNP
		if rec.is_snp and not rec.is_monomorphic:
			#Check if parsimony-informative
			if is_PIS(rec):
				# g=list()
				# for call in rec.samples:
				# 	if call.gt_bases:
				# 		g.append(call.gt_bases)
				# print(g)
				count+=1
				#if this is the final PIS, submit region to list
				if count == params.force:
					stop = rec.POS
					t = tuple([this_chrom, start, stop])
					regions.append(t)
					start = stop + 1
					count = 0

	t = tuple([this_chrom, start, contigs[this_chrom]])
	regions.append(t)

	print("Writing regions to out.regions...")
	write_regions("out.regions", regions)

#Function to write list of regions tuples, in GATK format
def write_regions(f, r):

	with open(f, 'w') as fh:
		try:
			for reg in r:
				ol = str(reg[0]) + ":" + str(reg[1]) + "-" + str(reg[2]) + "\n"
				fh.write(ol)
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()

#Function to check pyVCF record for if parsimony informative or not
def is_PIS(r):
	counts = dict()
	for sample in r.samples:
		if sample.gt_bases:
			nucs = [x.upper() for x in sample.gt_bases.split("/")]
			if nucs[0] == nucs[1] and nucs[0] in ["A", "C", "G", "T"]:
				if not nucs[0] in counts:
					counts[nucs[0]] = 1
				else:
					counts[nucs[0]] += 1
			else:
				for a in nucs:
					if a in ["A", "G", "C", "T"]:
						if not a in counts:
							counts[a] = 1
						else:
							counts[a] += 1
	has_two = 0
	if len(counts.keys()) > 1:
		for key in counts.keys():
			if counts[key] >= 2:
				has_two += 1
			if has_two >= 2:
				return(True)
	return(False)



#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'v:f:h', \
			["vcf=" "help", "force="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.force=100000

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
			elif opt in ('f','force'):
				self.force=int(arg)
			elif opt in ('h', 'help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Must provide VCF file <-v,--vcf>")

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfindBreaksVCF.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-v <input.vcf> -f <100000>\n")
		print ("Description: Breaks chromosomes into chunks of X parsimony-informative sites, for running MDL")

		print("""
	Arguments:
		-v,--vcf	: VCF file for parsing
		-f,--force	: Number of PIS to force a break
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
