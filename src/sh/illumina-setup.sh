#!/bin/bash
#
# Author: Peter Krusche <pkrusche@illumina.com>

# setup cmake environment vars for development on ilmn cluster
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ -d /illumina ]]; then
    if [[ "$(grep '6\.5' /etc/centos-release)" ]]; then
        . /etc/profile.d/modules.sh
    fi

    is_lua_modules=$(module --version 2>&1 | grep Lua)
    if [[ -z $is_lua_modules ]]; then
        unset MODULEPATH
        module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all
        module load CMake/3.2.1-GCC-4.9.2
        module load Java/1.8.0_40
    else
        ml purge  &> /dev/null
        ml HPCBIOS/2015q2
        ml sge/2011.11p1
        ml CMake/3.2.1-GCC-4.9.2
        ml Java/1.8.0_40
    fi
fi

export ANT_HOME=/illumina/thirdparty/ant/apache-ant-1.9.2
export EXTRA_CMAKE_OPTS="-DUSE_SGE=ON -DBUILD_VCFEVAL=ON -DVCFEVAL_WRAPPER=${DIR}/rtg-wrapper.sh"
export HAPPY_PATCH_SAMTOOLS="yes"
