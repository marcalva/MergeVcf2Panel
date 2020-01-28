#!/bin/bash

# Download htslib
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
tar xf htslib-1.10.2.tar.bz2
rm htslib-1.10.2.tar.bz2
cd htslib-1.10.2
./configure --prefix=$PWD
make
cd ../
