#!/bin/bash
# Wrapper for rtgtools to load Java first from a module

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ -d /illumina ]]; then
    if [[ "$(grep '6\.5' /etc/centos-release)" ]]; then
        . /etc/profile.d/modules.sh
    fi

    is_lua_modules=$(module --version 2>&1 | grep Lua)
    if [[ -z $is_lua_modules ]]; then
        unset MODULEPATH
        module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all
        module load Java/1.8.0_40
    else
        ml purge  &> /dev/null
        ml HPCBIOS/2015q2
        ml sge/2011.11p1
        ml Java/1.8.0_40
    fi
fi

${DIR}/rtg RTG_MEM=50g $@
