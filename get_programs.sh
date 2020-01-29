#!/bin/bash

# Download htslib 1.10.2
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
tar xf htslib-1.10.2.tar.bz2
rm htslib-1.10.2.tar.bz2
cd htslib-1.10.2
./configure --prefix=$PWD
make
cd ../

# Download bcftools
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
rm bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2
./configure --disable-lzma --disable-bz2 --prefix=$PWD
make
cd ../

