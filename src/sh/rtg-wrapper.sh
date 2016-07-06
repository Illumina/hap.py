#!/bin/bash
# Wrapper for rtgtools to load Java first from a module

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -e

module purge
unset MODULEPATH
module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all
module use /illumina/sync/software/thirdparty/HPCBIOS.2015q4/modules/all
module use ~/.local/easybuild/modules/all

module load Java/1.8.0_40

${DIR}/rtg RTG_MEM=50g $@
