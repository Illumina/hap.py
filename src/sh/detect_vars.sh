#!/bin/bash
#
# Detect common shell variables for testing
#
# Author: Peter Krusche <pkrusche@illumina.com>

if [[ -d "/illumina/thirdparty/graphviz/latest/bin" ]]; then
	export PATH=$PATH:/illumina/thirdparty/graphviz/latest/bin
fi

# fallback HG19 locations
if [[ ! -f "$HG19" ]]; then
	# HGREF variable can point to reference (hap.py asks for this, see src/python/Tools/__init__.py)
	export HG19="${HGREF}"
fi

if [[ ! -f "$HG19" ]]; then
	export HG19=/illumina/development/iSAAC/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
fi

if [[ ! -f "$HG19" ]]; then
	export HG19=~/workspace/human_genome/hg19.fa
fi

# fallback for finding hc binaries
if [[ ! -f "${HCDIR}/hap.py" ]]; then
	export HCDIR="$(pwd)/bin"
fi

if [[ ! -f "${HCDIR}/hap.py" ]]; then
	echo "Cannot find HC binaries. Set HCDIR to the bin directory."
	exit 1
fi

export PATH="$PATH:${HCDIR}"

# htslib is a shared object
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${HCDIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HCDIR}/lib

# detect or use python
if [[ -z ${PYTHON} ]]; then
	DEFAULT_PYTHON=1
	if [ -f "/illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh" ]; then
	    export PYTHON=/illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh
	fi
else
	DEFAULT_PYTHON=0
fi

export PYTHON=${PYTHON:-python}

PYVERSION=$(${PYTHON} --version 2>&1)
if [[ "$PYVERSION" != "Python 2.7."* ]] && [[ $DEFAULT_PYTHON == 1 ]]; then
	PYTHON=python2.7
fi

PYVERSION=$(${PYTHON} --version 2>&1)
if [[ "$PYVERSION" != "Python 2.7."* ]]; then
    echo "Hap.py requires Python 2.7.x. $PYTHON is $PYVERSION"
    exit 1
fi

export HCVERSION=`${PYTHON} ${HCDIR}/hap.py --version`
if [[ "$HCVERSION" == "" ]]; then
    echo "Cannot run hap.py to extract version information!"
    exit 1
fi
