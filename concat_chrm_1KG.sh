#!/bin/bash

# Cocnatenate chromosomes from 1KG into one file

bcftools="bcftools-1.10.2/bcftools"
tabix="htslib-1.10.2/tabix"

$bcftools concat \
    -O z \
    1KG/ALL.chr{{1..22},X}_GRCh38_sites.20170504.f.vcf.gz \
    > 1KG/ALL.GRCh38_sites.20170504.f.vcf.gz

$tabix -fp vcf 1KG/ALL.GRCh38_sites.20170504.f.vcf.gz

