#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export BOOST_ROOT=/illumina/thirdparty/boost/el6/boost_1_55_0_python2.6

# use GCC 4.8.2 + matching binutils and gdb
export PATH=/illumina/development/haplocompare/hc-virtualenv/bin:/illumina/thirdparty/gcc/el6/gcc-4.8.1/bin:/illumina/thirdparty/cmake/cmake-2.8.12.2/bin:$PATH

export LD_LIBRARY_PATH=/illumina/thirdparty/gcc/el6/gcc-4.8.1/lib

export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON"