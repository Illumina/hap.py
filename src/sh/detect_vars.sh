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
	export HCDIR="/illumina/development/haplocompare/haplocompare-devel/bin"
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
export PYTHON=${PYTHON:-python}
if [[ -z ${PYTHON} ]]; then
	if [ -d "/illumina/development/haplocompare/hc-virtualenv/bin" ]; then
	    export PYTHON=/illumina/development/haplocompare/hc-virtualenv/bin/python
	fi
fi

export HCVERSION=`${PYTHON} ${HCDIR}/hap.py --version`
