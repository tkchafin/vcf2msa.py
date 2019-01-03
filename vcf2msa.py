#!/usr/bin/python

import re
import sys
import subprocess
import os
import getopt
import vcf
import Bio
from os import path
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

def main():
	params = parseArgs()

	#Grab reference sequence first
	reference = dict()
	for contig in read_fasta(params.ref):
		if params.region:
			if contig[0].split()[0] == params.region.chr:
				name = contig[0].split()[0]
				fout = "contig_"+str(params.region.chr)+"_"+str(params.region.start)+"-"+str(params.region.end)+".fasta"
				print("Checking if",fout,"exists")
				if path.exists(fout) and params.force==False:
					print("Output file for contig already exists, skipping it:",fout)
				else:
					reference[name] = contig[1] #don't forget: 0-based index here, 1-based in VCF
		else:
			name = contig[0].split()[0]
			fout = "contig_"+str(name)+".fasta"
			print("Checking if",fout,"exists")
			if path.exists(fout) and params.force==False:
				print("Output file for contig already exists, skipping it:",fout)
			else:
				reference[name] = contig[1] #don't forget: 0-based index here, 1-based in VCF
	#read in samples
	vfh = vcf.Reader(filename=params.vcf)
	samples = list()
	sampleMask = dict() #dict of dict of sets
	for samp in vfh.samples:
		s = samp.split(".")[0]
		if s not in samples:
			samples.append(s)
		# if s not in sampleMask:
		# 	sampleMask[s] = dict()

	print("Found samples:",samples)

	#Get mask sites for each sample
	if len(reference) < 1:
		print("No contigs found.")
		sys.exit(0)

	for maskFile in params.mask:
		#print(maskFile)
		base=os.path.basename(maskFile)
		samp=base.split(".")[0]
		#print(samp)
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
	print("Found mask files:",sampleMask.keys())

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
				spos = 0
				epos = len(sequence)
				print("WARNING: Specified region outside of bounds: Using full contig",contig)
		#print(spos, " ", epos)

		#loop through each nuc position in contig
		for nuc in range(spos, epos):
			this_pos = dict()
			for samp in samples:
				this_pos[samp] = str()
			ref=None
			#For each sequence:
			#1. check if any samples have VCF data (write VARIANT)
			for rec in vfh.fetch(contig, nuc, nuc+1):
				#print(rec.samples)
				if int(rec.POS) != int(nuc+1):
					continue
					#sys.exit()
				if not ref:
					ref = rec.REF
				elif ref != rec.REF:
					print("Warning! Reference alleles don't match at position %s (contig %s): %s vs. %s" %(rec.POS, contig, ref, rec.REF))
				for ind in rec.samples:
					name = ind.sample.split(".")[0]
					if ind.gt_type:
						if not this_pos[name]:
							this_pos[name]=genotype_resolve(ind.gt_bases.split("/"),params.indel)
						else:
							this_pos[name]=genotype_resolve(ind.gt_bases.split("/"),params.indel, this_pos[name])

						#gt = "".join(sort(ind.gt_bases.split("/")))
						#print(ind.sample, " : ", ind.gt_bases)

			#3 insert Ns for masked samples at this position
			for samp in samples:
				if samp in sampleMask:
					if contig in sampleMask[samp] and nuc in sampleMask[samp][contig]:
						this_pos[samp] = "N"

			#4 if no allele chosen, write REF allele
			#use REF from VCF if possible, else pull from sequence
			for samp in samples:
				if not this_pos[samp] or this_pos[samp] == "":
					if not ref:
						this_pos[samp] = sequence[nuc]
					else:
						this_pos[samp] = ref

			#5. if there are insertions,  perform alignment to insert gaps
			l=None
			align=False
			for key in this_pos:
				if not l:
					l = len(this_pos[key])
				else:
					if l != len(this_pos[key]):
						#otherwise, set for alignment
						align=True

			#print(this_pos)
			if align==True:
				#print("aligning")
				#print(this_pos)
				this_pos = muscle_align(this_pos)
				#print(this_pos)

			#6 replace indels with Ns if they were masked
			maxlen = 1
			for key in this_pos:
				if len(this_pos[key]) > maxlen:
					maxlen = len(this_pos[key])
			if maxlen > 1:
				for samp in samples:
					if samp in sampleMask:
						if contig in sampleMask[samp] and nuc in sampleMask[samp][contig]:
							new = repeat_to_length("N", maxlen)
							this_pos[samp] = new
					elif samp not in this_pos:
						#sample was a multiple-nucleotide deletion
						new = repeat_to_length("-", maxlen)
						this_pos[samp] = new

			#7 Make sure nothing wonky happened
			p=False
			pis=False
			l=None
			a=None
			for key in this_pos:
				if not l:
					l = len(this_pos[key])
				else:
					if l != len(this_pos[key]):
						p=True
				# if not a:
				# 	if this_pos[key] != "N":
				# 		a = this_pos[key]
				# else:
				# 	if this_pos[key] != "N" and this_pos[key] != a:
				# 		pis = True
			if p==True:
				print("Warning: Alleles not of same length!!!")
				print(this_pos)
			# if pis==True:
			# 	print(this_pos)

			#8 add new bases to output string for contig/region
			for samp in samples:
				outputs[samp] += this_pos[samp]


		outFas = "contig_"+str(contig)+".fasta"
		if params.region:
			outFas = "contig_"+str(params.region.chr)+"_"+str(params.region.start)+"-"+str(params.region.end)+".fasta"
		with open(outFas, 'w') as fh:
			try:
				for sample in outputs:
					to_write = ">" + str(sample) + "\n" + outputs[sample] + "\n"
					fh.write(to_write)

			except IOError as e:
				print("Could not read file:",e)
				sys.exit(1)
			except Exception as e:
				print("Unexpected error:",e)
				sys.exit(1)
			finally:
				fh.close()


def repeat_to_length(string_to_expand, length):
    return (string_to_expand * (int(length/len(string_to_expand))+1))[:length]

#return dict alignment from dict of sequences, run via MUSCLE
def muscle_align(aln):
	records = Bio.Align.MultipleSeqAlignment([])
	for key, seq in aln.items():
		records.add_sequence(key, seq)

	cline = MuscleCommandline(clwstrict=True)
	child = subprocess.Popen(str(cline),
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		universal_newlines=True,
		shell=(sys.platform!="win32"))

	#write alignment to stdin handle for child process
	SeqIO.write(records, child.stdin, "fasta")
	child.stdin.close()

	#read alignment from stdout
	align = AlignIO.read(child.stdout, "clustal")

	new = dict()
	for record in align:
		new[record.id] = str(record.seq)
	return(new)



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

#internal function to resolve genotypes from VCF file
def genotype_resolve(l, indelPriority, e=None):
	if e:
		l.append(e)
	var=set()
	indel=set()
	if len(l) == 1:
		return(l[0])
	else:
		for gt in l:
			if len(gt) > 1 or gt=="*":
				if gt=="*":
					indel.add("-")
				else:
					indel.add(gt)
			else:
				var.add(gt)

		if indelPriority:
			if len(indel) == 1:
				return(next(iter(indel)))
			if len(indel) > 1:
				minlen=int()
				for i in indel:
					if not minlen:
						minlen = len(i)
					else:
						if minlen > len(i):
							minlen = len(i)
				for i in indel:
					if minlen == len(i):
						return(i)
		if len(var) == 1:
			return(next(iter(var)))
		elif len(var) > 1:
			gl = list(var)
			gl.sort()
			gt = "".join(gl)
			return(reverse_iupac_case(gt))
		elif len(var) < 1:
			if len(indel)==1:
				return(next(iter(indel)))
			elif len(indel) > 1:
				minlen=int()
				for i in indel:
					if not minlen:
						minlen = len(i)
					else:
						if minlen > len(i):
							minlen = len(i)
				for i in indel:
					if minlen == len(i):
						return(i)


#Function to translate a string of bases to an iupac ambiguity code, retains case
def reverse_iupac_case(char):
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N',
		'a':'a',
		'n':'n',
		'c':'c',
		'g':'g',
		't':'t',
		'ag':'r',
		'ct':'y',
		'ac':'m',
		'gt':'k',
		'at':'w',
		'cg':'s',
		'cgt':'b',
		'agt':'d',
		'act':'h',
		'acg':'v',
		'acgt':'n'
	}
	return iupac[char]


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
			["vcf=", "help", "ref=", "mask=","cov=","reg=", "indel", "force"])
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
		self.indel=False
		self.force=False

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
			if opt == "v" or opt == "vcf":
				self.vcf = arg
			elif opt == "c" or opt == "cov":
				self.cov=int(arg)
			elif opt == "m" or opt == "mask":
				self.mask.append(arg)
			elif opt =="r" or opt == "ref":
				self.ref=arg
			elif opt == "R" or opt == "reg":
				self.region = ChromRegion(arg)
			elif opt == "h" or opt == "help":
				pass
			elif opt == 'indel':
				self.indel=True
			elif opt == "force":
				self.force=True
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
		--indel		: In cases where indel conflicts with SNP call, give precedence to indel
		--force		: Overwrite existing alignment files
		-h,--help	: Displays help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
