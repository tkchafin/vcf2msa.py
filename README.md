# vcf2msa.py - Extracting multiple sequence alignments from VCF

Converts a multiple-sample VCF file to a multiple sequence alignment. Supports masking sites based on low coverage, by providing per-sample mpileup files. For these sites, we may not have sufficient information to assume that samples do not vary from the reference, and would like to instead write an N. Also supports indels, with indel realignment (to ascertain gap placement) performed in MUSCLE. 

This is meant to essentially perform the task of GATK's FastaAlternateReferenceMaker, but for multiple samples at once, in order to facilitate phylogenetic analysis. 

I've also included some associated scripts.

The code contained in this repository is provided for free via the GPL license, and is not guaranteed in any way. Use at your own risk.

### Dependencies
- Python >3
- pyVCF
- pySAM
- BioPython
- Muscle

The easiest way to install the dependencies is through conda:
```
conda install -c bioconda biopython pyvcf pysam muscle
```

### Inputs

vcf2msa.py requires two types of files: 1) A multi-sample VCF file containing your high-quality genotype calls (which may include indels); 2) an mpileup file for each sample, for every site (i.e. via running the samtools mpileup tool with the -aa option); and 3) a reference genome as FASTA (NOTE: this must be the referene used to build your VCF and mpileup files). 

You will additionally need to block-compress and index your VCF file:
```
#Run bgzip (NOT gzip) to compress your joint VCF file
bgzip file.vcf

#tabix to index it
tabix -h -f -p vcf file.vcf.gz
```

Note that you technically can run the script without the mpileup files, but I strongly warn against it because it means you are willing to assume that all samples share the REF allele, even if there is NO DATA (=no reads) to support that. Use at your own risk! :)

A note on sample names: In my GATK pipeline, I end up with a final 'joint variants' VCF file which contains two columns per sample: SampleID.variant and SampleID.variant2, with one containing the filtered SNP calls and the other the filtered indel calls. As a result, vcf2msa.py retains sample IDs as only the string preceeding the first "." and strips the remaining characters. So, for example if you have a sample named "s14A-B0.SNPs.filtered.calls", vcf2msa.py will only keep "s14A-B0" and treat any other samples with this prefix as identical. This might not be the desired behavior for you. This would be easy to change- if you need help with altering the code let me know and I can point you to the lines that need changing. 

Another consequence of how I treat sample IDs is that the per-sample mpileup files must ALSO be named in the same way, with the sample ID (i.e. everything preceding the first "." in the VCF column) should precede the first "." in the mpileup file. As a rule, I generally name these files as "sampleID.mpileup". This sampleID must be identical to that found in the VCF file, or else coverage information for that sample will not be incorporated. 

### Usage

View the help menu like so:
```
uaf90168:vcf2msa tkchafin$ python3 ./vcf2msa.py -h

Exiting because help menu was called.

vcf2msa.py

Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu
Description: Builds multiple sequence alignments from multi-sample VCF

	Mandatory arguments:
		-f,--fasta	: Reference genome file (FASTA)
		-v,--vcf	: VCF file containing genotypes

	Masking arguments:
		-c,--cov	: Minimum coverage to call a base as REF [default=1]
		-m,--mpileup: Per-sample mpileup file (one for each sample) for ALL sites (-aa in samtools)
		-d,--dp		: Toggle on to get depth from per-sample DP scores in VCF file

	Region selection arguments:
		Must use one of:
		-r,--reg	: Region to sample (e.g. chr1:1-1000)
		-R,--regfile: Text file containing regions to sample (1 per line, e.g. chr1:1-1000)
		-g,--gff	: Input GFF file containing regions to sample
			Note: All regions in the input file will be sampled, so this file
			will need to first be filtered to those regions you wish to select.
			Output files will be names according to the GFF name field

	Other arguments:
		-F,--flank	: Integer representing number of bases to add on either sides of selected regions [default=0]
		--id_field	: Field in GFF attributes (last column) giving the identifier to keep
				NOTE: Can list multiple (to be concatenated) like so: --id_field ID,Name
		--indel		: In cases where indel conflicts with SNP call, give precedence to indel
		--force		: Overwrite existing alignment files
		-h,--help	: Displays help menu
```
vcf2msa.py can either build an MSA file for each contig, or you can specify a region. I recommend specifying regions as parsing the whole files will take a long time. For example, if you wanted to fetch an MSA for chromosome 1 (chr1) positions 2050-10500, you could do:

```
#NOTE: No spaces allowed in the "-R" call
python3 ./vcf2msa.py -f <reference.fasta> -v <joint.vcf.gz> -r chr1:2050-10500
```

As I've mentioned, you should ideally specify per-sample mpileup files, and only call bases above some threshold coverage, and treat low (or null) coverage sites as ambiguous. To only call bases which have a coverage of 5 or greater, you could do:

```
#NOTE: No spaces allowed in the "-R" call
python3 ./vcf2msa.py -f <reference.fasta> -v <joint.vcf.gz> -r chr1:2050-10500 -m sample1.mpileup -m sample2.mpileup -m sample3.mpileup <...> -c 5
```
When doing this, you may find that vcf2msa.py runs very slow with large mpileup files- that is because I parse these files line by line to 1) find sites contained within the desired region; and 2) which fail the coverage threshold. To make this faster, you could use some bash commands to parse it down, or even do the coverage filtering yourself beforehand:
```
#just strip sites on wrong chromosome
grep "chr1" sample1.mpileup > sample1.subset.mpileup

#strip sites outside of the desired region:
grep "chr1" sample1.mpileup | awk '$2 > 2050 && $2 < 10500{print $0}' > sample1.subset.mpileup

#filter based on coverage
awk '$2 < 5{print $0}' sample1.mpileup > sample1.subset.mpileup
```

One final note: In cases where individuals are heterozygous for an indel AND a single-nucleotide substitution, the default behavior is to retain the SNP and ignore the indel. You can change this by setting an indel priority with <--indel>. Another default behavior is if an individual is heterozygous for two indels of different lengths, the shorter will always be retained. There is not currently a way built-in to change this.

In cases of >1 base indels, vcf2msa.py will locally re-align around the indel, to determine optimal gap placement, using MUSCLE. If this is the case, you will also need to have MUSCLE installed:

```
conda install -c bioconda muscle
```

# findBreaksVCF.py

You can use this script as a first-pass to delimit regions to extract alignments. In my pipeline, I find 'recombinatorial genes', or blocks of homogenous topology, using the MDL approach (see https://github.com/nstenz/TICR#mdl). This method finds breakpoints in a long (e.g. chromosomal) alignment using parsimony. To first create alignments for this, I recommend splitting your alignments into 'regions' representing a large number of parsimony-informative sites (since with the MDL approach these are the only sites carrying information). This is necessary to ease computation, since the number of sites with chromosomal and whole-genome alignments can easily get very large. 

This script defined 'breakpoints' in your VCF file as determined by a hard limit on number of parsimony-informative sites. 

```
uaf90168:vcf2msa tkchafin$ python3 ./findBreaksVCF.py 

Must provide VCF file <-v,--vcf>

findBreaksVCF.py

Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu

Usage:  ./findBreaksVCF.py -v <input.vcf> -f <100000>

Description: Breaks chromosomes into chunks of X parsimony-informative sites, for running MDL

	Arguments:
		-v,--vcf	: VCF file for parsing
		-f,--force	: Number of PIS to force a break
		-h,--help	: Displays help menu
```

To delimit breaks every 100,000 parsimony-informative sites, you would run it like so:

```
python3 ./findBreaksVCF.py -v joint.vcf -f 100000
```

The script will output a 'regions' file like so:
```
chr1:1-53445
chr1:53446-120054
...
...
...
```

