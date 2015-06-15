#!/bin/bash

# Simple performance and consistency test.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "GiaB / RTG test for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/GiaB/Complex_2ormoreindels_framerestoring_NIST2.19.ucsccoding.vcf \
			 	${DIR}/../../example/GiaB/Complex_2ormoreindels_framerestoring_RTG.ucsccoding.vcf \
			 	--fixchr-truth --fixchr-query \
			 	-o ${TMP_OUT} -P \
			 	-X \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

diff ${TMP_OUT}.counts.csv ${DIR}/../../example/GiaB/expected.counts.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.counts.csv ${DIR}/../../example/GiaB/expected.counts.csv"
	exit 1
fi

