#!/bin/bash

set -e

# Find python
PYTHON=python
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$1" == "clean" ]; then
    rm -rf zlib-1.2.8
    rm -rf boost_subset_1_58_0
    rm -rf boost_install
    rm -rf htslib
    rm -rf samtools
    rm -rf bcftools
    rm -rf rtg-tools
    exit 0
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

if [ -z "$BOOST_ROOT" ];
then
    if [ ! -d boost_install ];
    then
        tar xjf boost_subset_1_58_0.tar.bz2
        cd boost_subset_1_58_0
        ./bootstrap.sh
        ./b2 link=static -j4 --prefix=$DIR/boost_install -sZLIB_SOURCE=$DIR/zlib-1.2.8
        ./b2 link=static -j4 --prefix=$DIR/boost_install install -sZLIB_SOURCE=$DIR/zlib-1.2.8
        cd ..
    else
        echo "Boost already built. To rebuild, delete external/boost_install"
    fi
else
    echo "BOOST_ROOT is set, not building boost."
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
    # Samtools 1.3 seems to work out of the box
    # if [[ "$HAPPY_PATCH_SAMTOOLS" == "yes" ]]; then
    #     echo "Patching Samtools for compatibility with EasyBuild"
    #     patch -p0 samtools/Makefile < SAMtools_Makefile.patch
    # fi
    make -j4 -C samtools
else
    echo "samtools already built. To rebuild, delete external/samtools"
fi

# get vcfeval
# https://github.com/RealTimeGenomics/rtg-tools/archive/ga4gh-test.zip
if [[ ! -z $BUILD_VCFEVAL ]]; then
    if [[ ! -d rtg-tools/rtg-tools-install ]]; then
        echo "Building rtg-tools"
        mkdir -p rtg-tools
        cd rtg-tools
        wget http://github.com/RealTimeGenomics/rtg-tools/archive/296c61ed18e363574fdbc982bbe73c0b86c796ce.tar.gz -O rtg-tools.tar.gz
        tar xvf rtg-tools.tar.gz
        cd rtg-tools-296c61ed18e363574fdbc982bbe73c0b86c796ce

        if [[ ! -z ${ANT_HOME} ]]; then
            $ANT_HOME/bin/ant zip-nojre
        else
            ant zip-nojre
        fi
        cd ..

        RTG_ZIPFILE=$(ls rtg-tools-296c61ed18e363574fdbc982bbe73c0b86c796ce/dist/*-nojre.zip | head -1)
        RTG_BASE=$(basename $RTG_ZIPFILE -nojre.zip)
        jar xvf $RTG_ZIPFILE
        mv $RTG_BASE rtg-tools-install
        cp ../rtg.cfg rtg-tools-install
        chmod +x rtg-tools-install/rtg
        cd ..
    else
        echo "rtg-tools is already built. To rebuild, delete external/rtg-tools"
    fi
    if [[ -f ${VCFEVAL_WRAPPER} ]]; then
        cd rtg-tools
        echo "using wrapper for rtg-tools: ${VCFEVAL_WRAPPER}"
        cp ${VCFEVAL_WRAPPER} rtg-tools-install/rtg-wrapper.sh
        chmod +x rtg-tools-install/rtg
    fi
fi
