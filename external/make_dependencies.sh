#!/bin/bash

set -e

# Find python
PYTHON=python
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -d "/illumina/thirdparty/Python-2.7.3/bin" ]; then
    PYTHON=/illumina/thirdparty/Python-2.7.3/bin/python
fi
if [ -d "/illumina/development/pgtools/pgt-virtualenv/bin" ]; then
    PYTHON=/illumina/development/pgtools/pgt-virtualenv/bin/python
fi

if [ ! -f zlib-1.2.8/libz.a ];
then
    tar xzf zlib-1.2.8.tar.gz
    cd zlib-1.2.8
    ./configure
    cd ..
    make -j4 -C zlib-1.2.8
else
    echo "Zlib already built. To rebuild, delete external/zlib-1.2.8"
fi

# Get from git
#
# git clone git@github.com:samtools/htslib.git
# git clone git@github.com:samtools/bcftools.git
# git clone git@github.com:samtools/samtools.git

if [ ! -f htslib/libhts.so ] || [ -f htslib/libhts.dyld ] ;
then
    tar xzf htslib.tar.gz
    make -j4 -C htslib

    # windows shared folder workaround
    if [ -e "htslib/libhts.so" ] && [ ! -e "htslib/libhts.so.1" ] ;
    then
        cp htslib/libhts.so htslib/libhts.so.1
    fi
else
    echo "HTSLIB already built. To rebuild, delete external/htslib"
fi

if [ ! -f bcftools/bcftools ];
then
    tar xzf bcftools.tar.gz
    make -j4 -C bcftools
else
    echo "bcftools already built. To rebuild, delete external/bcftools"
fi

if [ ! -f samtools/samtools ];
then
    tar xzf samtools.tar.gz
    make -j4 -C samtools
else
    echo "samtools already built. To rebuild, delete external/samtools"
fi

if [[ "$1" == "build_muscle" ]]; then
    if [ ! -f muscle3.8.31/src/muscle ];
    then
        wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz
        tar xzf muscle3.8.31_src.tar.gz
        make -j4 -C muscle3.8.31/src
    else
        echo "muscle already built. To rebuild, delete external/muscle3.8.31"
    fi
else
    echo "Muscle3 is disabled. Set the USE_MUSCLE variable to enable."
fi

