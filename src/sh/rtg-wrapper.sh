#!/bin/bash
# Wrapper for rtgtools to load Java first from a module

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -e

if [[ -d /illumina ]]; then
    if [[ "$(grep '6\.5' /etc/centos-release)" ]]; then
        . /etc/profile.d/modules.sh
    fi

    module purge --force &> /dev/null

    is_lua_modules=$(module --version 2>&1 | grep Lua)
    if [[ -z $is_lua_modules ]]; then
        unset MODULEPATH
        module use /illumina/sync/software/unofficial/HPCBIOS/2016q2/modules/all &> /dev/null
    else
        module purge  &> /dev/null
        module load --force HPCBIOS/2016q2 &> /dev/null
        module load sge/2011.11p1
    fi

    module load Java/1.8.0_40
fi

${DIR}/rtg RTG_MEM=50g $@
