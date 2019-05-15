#!/usr/bin/python

import re
import sys
import os
import getopt

def main():
	params = parseArgs()




#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 's:e:f:hz', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.start=None
		self.int=None
		self.zero=False

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
		print ("\nregionsFromFasta.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-f /path/to/fasta -s <int> -e <int>\n")
		print ("Description: Extracts a subalignment from a multi-fasta, given start and stop coordiantes")

		print("""
	Arguments:
		-f	: Fasta file
		-s	: Start coordinate (default: 1st base)
		-e	: End coordinate (default: last base)
		-z	: (Boolean) toggle on if coordiantes are zero-based indexing
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
