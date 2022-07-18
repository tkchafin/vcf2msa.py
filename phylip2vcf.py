#!/usr/bin/python

import re
import sys
import subprocess
import os
import getopt
from collections import Counter


def main():

	header="##fileformat=VCFv4\n"
	header=header + "##Spoofed VCF created by phylip2vcf.py\n"
	header=header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

	contig="SPOOF"
	id="."
	qual="999"
	filter="."
	info="."
	format="GT"

	params = parseArgs()

	dat = read_phylip(params.phy)

	# add samples to header
	for k in dat.keys():
		header = header + "\t" + str(k)
	header = header + "\n"
	#print(header)

	records=list()
	for i in range(0, len(dat[list(dat.keys())[0]])):
		rec=contig + "\t" + str(i+1) + "\t" + id + "\t"
		alleles=list()
		for ind in dat.keys():
			alleles.extend(get_iupac_caseless(dat[ind][i]))
		c=Counter(alleles)
		oc=c.most_common()
		ref=None
		alts=list()
		map=dict()
		index=0
		map["N"]="."
		for n in oc:
			if n[0] != "N":
				if not ref:
					ref=n[0]
					map[ref]=str(index)
					index+=1
				else:
					alts.append(n[0])
					map[n[0]]=str(index)
					index+=1
		# if data was all missing, set fake ref
		if not ref:
			ref="A"
			map[ref]="0"
		# if data was monomorphic, set fake alt
		if len(alts) == 0:
			alts="T"
			map["T"] = "1"

		rec = rec + ref + "\t" + str(",".join(alts)) + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format
		#print(rec)

		for k in dat.keys():
			gt=get_vcf_genotype(dat[k][i], map)
			rec = rec + "\t" + gt
		rec = rec + "\n"
		records.append(rec)
		#print(rec)
		#sys.exit()

	with open(params.vcf, "w") as ofh:
		ofh.write(header)
		for rec in records:
			ofh.write(rec)


def get_vcf_genotype(gt, map):
	alleles=get_iupac_caseless(gt)
	ret=map[alleles[0]] + "/" + map[alleles[1]]
	return(ret)

def get_major_allele(l, num=None, vcf=False):
    """Get most common alleles in list.

    Args:
        l (List[str]): List of genotypes for one sample.

        num (int, optional): Number of elements to return. Defaults to None.

        vcf (bool, optional): Alleles in VCF or STRUCTURE-style format. Defaults to False.

    Returns:
        list: Most common alleles in descending order.
    """
    all_items = list()
    for i in l:
        if vcf:
            all_items.extend(i.split("/"))
        else:
            all_items.extend(get_iupac_caseless(i))

    c = Counter(all_items)  # requires collections import
    rets = c.most_common(num)

    # Returns two most common non-ambiguous bases
    # Makes sure the least common base isn't N or -9
    if vcf:
        return [x[0] for x in rets if x[0] != "-9"]
    else:
        return [x[0] for x in rets if x[0] in ["A", "T", "G", "C"]]


def count_alleles(l, vcf=False):
    """Count how many total alleles there are.

    Args:
        l (List[str]): List of IUPAC or VCF-style (e.g. 0/1) genotypes.
        vcf (bool, optional): Whether genotypes are VCF or STRUCTURE-style. Defaults to False.

    Returns:
        int: Total number of alleles in l.
    """
    all_items = list()
    for i in l:
        if vcf:
            all_items.extend(i.split("/"))
        else:
            all_items.extend(get_iupac_caseless(i))
    all_items = remove_items(all_items, ["-9", "-", "N", -9])
    return len(set(all_items))

def read_phylip(phy):
	if os.path.exists(phy):
		with open(phy, 'r') as fh:
			try:
				num=0
				ret = dict()
				for line in fh:
					line = line.strip()
					if not line:
						continue
					num += 1
					if num == 1:
						continue
					arr = line.split()
					ret[arr[0]] = list(arr[1])
				return(ret)
			except IOError:
				print("Could not read file ",phy)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%phy)

def get_iupac_caseless(char):
    """Split IUPAC code to two primary characters, assuming diploidy.

    Gives all non-valid ambiguities as N.

    Args:
        char (str): Base to expand into diploid list.

    Returns:
        List[str]: List of the two expanded alleles.
    """
    lower = False
    if char.islower():
        lower = True
        char = char.upper()
    iupac = {
        "A": ["A", "A"],
        "G": ["G", "G"],
        "C": ["C", "C"],
        "T": ["T", "T"],
        "N": ["N", "N"],
        "-": ["N", "N"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["N", "N"],
        "D": ["N", "N"],
        "H": ["N", "N"],
        "V": ["N", "N"],
    }
    ret = iupac[char]
    if lower:
        ret = [c.lower() for c in ret]
    return ret

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


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hp:v:', \
			["vcf=", "help", "phylip=", "phy="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.phy=None

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
			elif opt == "p" or opt=="phy" or opt=="phylip":
				self.phy=arg
			elif opt == "h" or opt == "help":
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Must provide VCF output name <-v,--vcf>")
		if not self.phy:
			self.display_help("Must provide phylip input file <-p,--phy>")

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nphylip2vcf.py\n")
		print ("Contact:Tyler K. Chafin, tylerkchafin@gmail.com")
		print ("Description: Converts phylip to spoof VCF file")

		print("""
	Mandatory arguments:
		-p,--phylip	: Phylip input file
		-v,--vcf	: VCF file containing genotypes

	Other arguments:
		-h,--help	: Displays help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
