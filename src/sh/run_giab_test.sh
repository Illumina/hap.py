#!/bin/bash

# Simple performance and consistency test.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Small GiaB / RTG test for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/GiaB/Complex_2ormoreindels_framerestoring_NIST2.19.ucsccoding.vcf \
			 	${DIR}/../../example/GiaB/Complex_2ormoreindels_framerestoring_RTG.ucsccoding.vcf \
			 	-o ${TMP_OUT} \
			 	-X \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

echo "Large GiaB / RTG test for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/NIST_indels/Complex_1ormoreindels_NIST2.19.vcf.gz \
			 	${DIR}/../../example/NIST_indels/Complex_1ormoreindels_RTG.vcf.gz \
			 	-o ${TMP_OUT} -l chr21 \
			 	-X --verbose \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/NIST_indels/expected.summary.21.csv
if [[ $? != 0 ]]; then
	echo "summary differs! -- diff ${TMP_OUT}.summary.csv ${DIR}/../../example/NIST_indels/expected.summary.21.csv"
	exit 1
fi

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/NIST_indels/Complex_1ormoreindels_NIST2.19.vcf.gz \
			 	${DIR}/../../example/NIST_indels/Complex_1ormoreindels_RTG.vcf.gz \
			 	-o ${TMP_OUT} -l chr1 \
			 	-X --verbose \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/NIST_indels/expected.summary.1.csv
if [[ $? != 0 ]]; then
	echo "summary differs! -- diff ${TMP_OUT}.summary.csv ${DIR}/../../example/NIST_indels/expected.summary.1.csv"
	exit 1
fi

