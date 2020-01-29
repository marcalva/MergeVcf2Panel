#!/bin/bash

# Download variant sites from 1KG

bcftools="../bcftools-1.10.2/bcftools"
bgzip="../htslib-1.10.2/bgzip"
tabix="../htslib-1.10.2/tabix"

mkdir -p 1KG
cd 1KG
for i in $( seq 1 22 ) X Y; do
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${i}_GRCh38_sites.20170504.vcf.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${i}_GRCh38_sites.20170504.vcf.gz.tbi
    $bcftools view ALL.chr${i}_GRCh38_sites.20170504.vcf.gz | \
        cut -f 1-5 | \
        $bgzip -c > ALL.chr${i}_GRCh38_sites.20170504.f.vcf.gz
    $tabix -fp vcf ALL.chr${i}_GRCh38_sites.20170504.f.vcf.gz
done

touch *tbi

cd ../
