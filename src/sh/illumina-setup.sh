#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. /etc/profile.d/modules.sh
module purge
unset MODULEPATH
module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all
module load CMake/3.2.1-GCC-4.9.2

# for rtgtools / vcfeval
module load Java/1.8.0_40
export ANT_HOME=/illumina/thirdparty/ant/apache-ant-1.9.2

export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON -DBUILD_VCFEVAL=ON -DVCFEVAL_WRAPPER=${DIR}/rtg-wrapper.sh"
export HAPPY_PATCH_SAMTOOLS="yes"

