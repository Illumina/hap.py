#!/bin/bash
#
# Setup script for custom locations of Boost and Python
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ukch-dev
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ -d /illumina ]]; then
    if [[ "$(grep '6\.5' /etc/centos-release)" ]]; then
        . /etc/profile.d/modules.sh
    fi

    is_lua_modules=$(module --version 2>&1 | grep Lua)
    if [[ -z $is_lua_modules ]]; then
        unset MODULEPATH
        module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all &> /dev/null
        module use /illumina/sync/software/unofficial/HPCBIOS/2016q2/modules/all &> /dev/null
    else
        module purge  &> /dev/null
        module load --force HPCBIOS/2016q2 &> /dev/null
        module load sge/2011.11p1
    fi
fi

module load CMake/3.2.1-GCC-4.9.2
module load Java/1.8.0_40

export ANT_HOME=/illumina/thirdparty/ant/apache-ant-1.9.2
export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON -DBUILD_VCFEVAL=ON -DVCFEVAL_WRAPPER=${DIR}/rtg-wrapper.sh"
export HAPPY_PATCH_SAMTOOLS="yes"

