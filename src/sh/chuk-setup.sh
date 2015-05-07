#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export BOOST_ROOT=/illumina/thirdparty/boost/boost_1_55_0_python2.4

# use GCC 4.8.2 + matching binutils and gdb
export PATH=/illumina/development/haplocompare/hc-virtualenv/bin:/illumina/thirdparty/gcc/el5/gcc-4.9.2/bin:/illumina/thirdparty/gdb/gdb-7.8/bin:/illumina/thirdparty/binutils/binutils-2.24/bin:/illumina/thirdparty/cmake/el5/cmake-2.8.12.1/bin:$PATH

export PYTHONHOME=/illumina/thirdparty/python/python-2.7.5
unset PYTHONUSERBASE

export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON"