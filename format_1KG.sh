#!/bin/bash

# Merge the downloaded 1KG chromosome files into one file

vcfout="1KG/ALL.GRCh38_sites.20170504.vcf"

zcat 1KG/ALL.chr1_GRCh38_sites.20170504.vcf.gz | \
    cut -f 1-5 > $vcfout

for i in $( seq 2 22 ) X Y; do
    zcat 1KG/ALL.chr${i}_GRCh38_sites.20170504.vcf.gz | \
        awk 'substr($1, 1, 1) != "#"' | \
        cut -f 1-5 >> $vcfout
done

htslib-1.10.2/bgzip -f $vcfout

