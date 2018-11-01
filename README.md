# vcf2msa

Converts a multiple-sample VCF file to a multiple sequence alignment. Supports regions and masking sites based on coverage. For example, a VCF or mpileup file can be given which specifies sites with low coverage. For these sites, we may not have sufficient information to assume that samples do not vary from the reference, and would like to instead write an N. This is meant to essentially perform the task of GATK's FastaAlternateReferenceMaker, but for multiple samples at once, in order to facilitate phylogenetic analysis. 

The code contained in this repository is provided for free via the GPL license, and is not guaranteed in any way. Use at your own risk.

### Dependencies
- Python >3
- pyVCF
- pySAM
- BioPython

The easiest way to install the dependencies is through conda:
```
conda install -c bioconda biopython pyvcf pysam
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
		-h,--help	: Displays help menu
    
```
vcf2msa.py can either build an MSA file for each contig, or you can specify a region. I recommend specifying regions as parsing the whole files will take a long time. For example, if you wanted to fetch an MSA for chromosome 1 (chr1) positions 2050-10500, you could do:

```
#NOTE: No spaces allowed in the "-R" call
python3 ./vcf2msa.py -r <reference.fasta> -v <joint.vcf.gz> -R chr1:2050-10500
```

As I've mentioned, you should ideally specify per-sample mpileup files, and only call bases above some threshold coverage, and treat low (or null) coverage sites as ambiguous. To only call bases which have a coverage of 5 or greater, you could do:

```
#NOTE: No spaces allowed in the "-R" call
python3 ./vcf2msa.py -r <reference.fasta> -v <joint.vcf.gz> -R chr1:2050-10500 -m sample1.mpileup -m sample2.mpileup -m sample3.mpileup <...> -c 5
```

One final note: In cases where individuals are heterozygous for an indel AND a single-nucleotide substitution, the default behavior is to retain the SNP and ignore the indel. You can change this by setting an indel priority with <--indel>. Another default behavior is if an individual is heterozygous for two indels of different lengths, the shorter will always be retained. There is not currently a way built-in to change this.


# altRefMaker.py

This script is a naive replication of the GATK FastaAlternateReferenceMaker tool, and is in progress. Use at your own risk.


# To do:
- Figure out a way to index the mpileup files? Taking too long to parse
- Documentation for altRefMaker.py
