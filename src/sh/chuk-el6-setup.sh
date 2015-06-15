#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export BOOST_ROOT=/illumina/thirdparty/boost/boost_1_57_0_python2.7-fPIC

module purge
unset MODULEPATH
module use /illumina/sync/software/thirdparty/HPCBIOS.20150417/modules/all
module load CMake
module unload GCC
module load GCC/4.9.2
module unload zlib
module unload ncurses

export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON"
