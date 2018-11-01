# vcf2msa

Converts a multiple-sample VCF file to a multiple sequence alignment. Supports regions and masking sites. For example, a VCF or mpileup file can be given which specifies sites with low coverage. For these sites, we may not have sufficient information to assume that samples do not vary from the reference, and would like to instead write an N. 

In progress- not yet ready for human consumption. 

Documentation to come later :)

Requires:
- Python >3
- pyVCF
- pySAM

Need to BGZIP VCF and TABI-index it:
```
bgzip file.vcf
tabix -h -f -p vcf file.vcf.gz
```

Note: figure out a way to index the mpileup files? Taking too long to parse
