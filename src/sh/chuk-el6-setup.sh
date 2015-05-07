#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export BOOST_ROOT=/illumina/thirdparty/boost/boost_1_57_0_python2.7-fPIC

# use GCC 4.8.2 + matching binutils and gdb
export PATH=/illumina/development/haplocompare/hc-virtualenv/bin:/illumina/thirdparty/gcc/gcc-4.9.2/bin:/illumina/thirdparty/gdb/gdb-7.9/bin:/illumina/thirdparty/binutils/binutils-2.25/bin:/illumina/thirdparty/cmake/cmake-3.1.3/bin:$PATH

export PYTHONHOME=/illumina/thirdparty/python/python-2.7.8
unset PYTHONUSERBASE

export LD_RUN_PATH=$LD_LIBRARY_PATH:/illumina/thirdparty/gcc/gcc-4.9.2/lib64:/illumina/thirdparty/gcc/gcc-4.9.2/lib

export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON"
