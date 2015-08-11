#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

module purge
unset MODULEPATH
module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all
module load CMake/3.2.1-GCC-4.9.2
# These break samtools. We need to make a patch for this.

export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON"
export HAPPY_PATCH_SAMTOOLS="yes"
