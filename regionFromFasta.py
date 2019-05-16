#!/usr/bin/python

import re
import sys
import os
import getopt

def main():
	params = parseArgs()

	s=0
	e=-1
	if params.start:
		if params.zero == False:
			s = params.start
		else:
			if params.start <= 0:
				print("Start coordinate cannot be less than 1 unless zero-based indexing <-z>")
				sys.exit(1)
			s=params.start-1
	#print("start=",s)

	if params.end:
		if params.zero == False:
			e = params.end
		else:
			e=params.end-1
	#print("end=",e)
	newSamps = list()
	for seq in read_fasta(params.fasta):
		keep=None
		if e > len(seq[1]):
			print("Warning: end coordinate",e,"exceeds sequence length.")
			keep = seq[1][s:]
		else:
			keep=seq[1][s:e]
		if keep:
			newSamps.append((seq[0], keep))
		else:
			print("No sequence was kept. Something went wrong.")

	write_fasta(params.out,newSamps)

#Function to write fasta-formatted sequences
def write_fasta(f, aln):

	with open(f, 'w') as fh:
		try:
			for samp in aln:
				ol = ">" + str(samp[0]) + "\n" + str(samp[1]) + "\n"
				fh.write(ol)
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()

#Read samples as FASTA. Generator function
def read_fasta(fas):

	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
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
						split_line = line.split()
						contig = (split_line[0].replace(">",""))
					else:
						seq += line
				#Iyield last sequence, if it has both a header and sequence
				if contig and seq:
					yield([contig,seq])
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 's:e:f:hzo:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.start=None
		self.end=None
		self.zero=False
		self.fasta=None
		self.out="out.fasta"

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
			if opt == "s":
				self.start = int(arg)
			elif opt == "e":
				self.end = int(arg)
			elif opt == "z":
				self.zero=True
			elif opt == "f":
				self.fasta = arg
			elif opt == "o":
				self.out = arg
			elif opt in ('h', 'help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasta:
			self.display_help("Must provide FASTA file <-f>")

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nregionsFromFasta.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-f /path/to/fasta -s <int> -e <int>\n")
		print ("Description: Extracts a subalignment from a multi-fasta, given start and stop coordiantes")

		print("""
	Arguments:
		-f	: Fasta file
		-s	: Start coordinate (default: 1st base)
		-e	: End coordinate (default: last base)
		-z	: (Boolean) toggle on if coordiantes are zero-based indexing (default False)
		-o	: Output file name (default=out.fasta)
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
