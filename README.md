# genomeAlignMaker

Contains scripts to make nucleotide alignments from genomic polymorphism data. 

### altRefMaker.py
Pulls an 'alternate' reference sequence for all contigs given a VCF file and mpileup file representing depth of coverage for all bases. This is a naive re-implementation of the GATK FastaAlternateReferenceMaker tool.

### alignMaker.pt
Uses a sliding window down each contig to find regions of the genome which have suitable per-base and per-sample coverage for phylogenetic analysis.
