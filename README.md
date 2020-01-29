# Merge VCF file into a reference

This python script is designed to merge VCF files with that of a 
reference. This is useful for pre-processing multiple VCF files 
that may have been genotyped/sequenced on different platforms. This 
can result in distinct annotation of variants, for example the 
strand may be flipped.

Merging the VCF files involves:
* removing variants that are not present in the reference
* removing SNVs that are strand ambiguous
* flipping the strand
* flipping the reference and alternate alleles

The stand and alleles must be unambiguous in order for them to 
be flipped. For example, an A/T SNP is ambiguous. A/T can 
specify reference/alternate or forward/reverse strand.

## Usage

The usage is as follows:

```bash
python MergeVcf2Panel.py \
           --ref 1KG/ALL.chr1_GRCh38_sites.20170504.vcf.gz \
           --vcfin my_geno.vcf \
           --vcfout my_geno.chr1.1KG.vcf
```

This takes the reference file from `1KG/ALL.chr1_GRCh38_sites.20170504.vcf.gz`
, the first chromosome, and merges the variants from `my_geno.vcf` into 
`my_geno.chr1.1KG.vcf`.

## Reference

The reference file is a VCF where only the first 5 columns are used 
(corresponding to location, ID, and alleles). 
A popular reference panel is the 1,000 Genomes. The `download_1KG.sh` 
script downloads the variant sites in VCF format into the folder 
`1KG`. These are variants lifted over into GRCh38 coordinates. 
These files can be used as the reference input for the 
`MergeVcf2Panel.py` script. To merge the chromosomal files and 
remove the columns after the fifth, run `get_htslib.sh` to get 
the bgip tools and then run `format_1KG.sh`. This places all 
the variants into `1KG/ALL.GRCh38_sites.20170504.vcf`, which 
can be used with the python script.

If you want to run all chromosomes in a single serial run, merge 
VCF with `concat_chrm_1KG.sh`.
