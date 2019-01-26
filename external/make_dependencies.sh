#!/bin/bash

set -e

# Find python
PYTHON=python
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TLD=$(pwd)/scratch
ISD=$(pwd)

if [ "$1" == "rebuild" ];
then
    rm -rf ${TLD}
    rm -rf ${ISD}/include/boost
    rm -rf ${ISD}/include/htslib
    rm -rf ${ISD}/bin/bcftools
    rm -rf ${ISD}/bin/samtools
    rm -rf ${ISD}/libexec/rtg-tools-install
fi

mkdir -p ${TLD}
if [ ! -f ${TLD}/zlib-1.2.8/libz.a ];
then
    cd ${TLD}
    rm -rf ${TLD}/zlib-1.2.8
    tar xzf ${DIR}/zlib-1.2.8.tar.gz
    cd zlib-1.2.8
    ./configure --prefix ${ISD}
    make -j4
    make -j4 install
else
    echo "Zlib already built. To rebuild, delete ${TLD}/zlib-1.2.8"
fi

if [ -z "$BOOST_ROOT" ];
then
    if [ ! -d ${ISD}/include/boost ];
    then
        cd ${TLD}
        rm -rf ${TLD}/boost_subset_1_58_0
        tar xjf ${DIR}/boost_subset_1_58_0.tar.bz2
        cd boost_subset_1_58_0
        ./bootstrap.sh
        ./b2 link=static -j4 --prefix=$ISD -sZLIB_SOURCE=$TLD
        ./b2 link=static -j4 --prefix=$ISD install -sZLIB_SOURCE=$TLD/zlib-1.2.8
    else
        echo "Boost already built. To rebuild, delete ${ISD}/include/boost"
    fi
else
    echo "BOOST_ROOT is set, not building boost."
fi

if [ ! -d ${ISD}/include/htslib ] ;
then
    cd ${TLD}
    rm -rf ${TLD}/htslib
    tar xzf ${DIR}/htslib.tar.gz
    cd htslib
    ./configure --prefix=${ISD} \
        CFLAGS=-I${ISD}/include\ -g \
        CXXFLAGS=-I${ISD}/include\ -g \
        LDFLAGS=-L${ISD}/lib \
        --disable-plugins \
        --disable-libcurl \
        --disable-lzma \
        --disable-bz2
    make -j4
    make -j4 install

    # windows shared folder workaround
    if [ -e "${ISD}/lib/libhts.so" ] && [ ! -e "${ISD}/lib/libhts.so.1" ] ;
    then
        cp ${ISD}/lib/libhts.so ${ISD}/lib/libhts.so.1
    fi
else
    echo "HTSLIB already built. To rebuild, delete ${ISD}/include/htslib"
fi

if [ ! -f ${ISD}/bin/bcftools ];
then
    cd ${TLD}
    rm -rf ${TLD}/bcftools
    tar xzf ${DIR}/bcftools.tar.gz
    cd bcftools
    make -j4 prefix=${ISD}
    make -j4 prefix=${ISD} install
else
    echo "bcftools already built. To rebuild, delete ${ISD}/bin/bcftools"
fi

if [ ! -f ${ISD}/bin/samtools ];
then
    cd ${TLD}
    rm -rf ${TLD}/samtools
    tar xzf ${DIR}/samtools.tar.gz
    cd samtools
    autoconf -Wno-syntax || autoconf -Wno-syntax
    ./configure --prefix=${ISD} \
        --with-htslib=${TLD}/htslib \
        --without-curses \
        CFLAGS=-I${ISD}/include \
        CPPFLAGS=-I${ISD}/include \
        LDFLAGS=-L${ISD}/lib
    make -j4
    make -j4 install
else
    echo "samtools already built. To rebuild, delete ${ISD}/bin/samtools"
fi

# get vcfeval
# https://github.com/RealTimeGenomics/rtg-tools/archive/ga4gh-test.zip
if [[ ! -z $BUILD_VCFEVAL ]]; then
    if [[ ! -d ${ISD}/libexec/rtg-tools-install ]]; then
        echo "Building rtg-tools"
        cd ${TLD}
        mkdir -p ${TLD}/rtg-tools
        cd rtg-tools
        wget http://github.com/RealTimeGenomics/rtg-tools/archive/3.10.1.tar.gz -O ${TLD}/rtg-tools/rtg-tools.tar.gz
        tar xvf rtg-tools.tar.gz
        cd rtg-tools-3.10.1

        if [[ ! -z ${ANT_HOME} ]]; then
            $ANT_HOME/bin/ant zip-nojre
        else
            ant zip-nojre
        fi
        cd ..

        RTG_ZIPFILE=$(ls rtg-tools-3.10.1/dist/*-nojre.zip | head -1)
        RTG_BASE=$(basename $RTG_ZIPFILE -nojre.zip)
        jar xvf $RTG_ZIPFILE
        mv $RTG_BASE ${ISD}/libexec/rtg-tools-install
        cp ${DIR}/rtg.cfg ${ISD}/libexec/rtg-tools-install
        chmod +x ${ISD}/libexec/rtg-tools-install/rtg
    else
        echo "rtg-tools is already built. To rebuild, delete ${ISD}/libexec/rtg-tools-install"
    fi
    if [[ -f ${VCFEVAL_WRAPPER} ]]; then
        echo "using wrapper for rtg-tools: ${VCFEVAL_WRAPPER}"
        cp ${VCFEVAL_WRAPPER} ${ISD}/libexec/rtg-tools-install/rtg-wrapper.sh
        chmod +x ${ISD}/libexec/rtg-tools-install/rtg-wrapper.sh
    fi
fi
