#!/bin/bash

# Test if chr prefixes are detected correctly automatically
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Test for chr prefix detection  ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/numeric_chrs/truth.vcf \
			 	${DIR}/../data/numeric_chrs/query.vcf \
                -f ${DIR}/../data/numeric_chrs/fp.bed \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/numeric_chrs/test.fa \
                -V \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../data/numeric_chrs/expected.summary.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.summary.csv ${DIR}/../data/numeric_chrs/expected.summary.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz | grep -v ^# > ${TMP_OUT}.vcf
diff ${TMP_OUT}.vcf ${DIR}/../data/numeric_chrs/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../data/numeric_chrs/expected.vcf"
	exit 1
else
    echo "Numeric chr test successful"
    rm -rf ${TMP_OUT}.*
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/numeric_chrs/chrtruth.vcf \
			 	${DIR}/../data/numeric_chrs/chrquery.vcf \
                -f ${DIR}/../data/numeric_chrs/chrfp.bed \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/numeric_chrs/chrtest.fa \
                -V \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../data/numeric_chrs/expected.summary.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.summary.csv ${DIR}/../data/numeric_chrs/expected.summary.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz | grep -v ^# > ${TMP_OUT}.vcf
diff ${TMP_OUT}.vcf ${DIR}/../data/numeric_chrs/chrexpected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../data/numeric_chrs/chrexpected.vcf"
	exit 1
else
    echo "Chr prefixed test successful"
    rm -rf ${TMP_OUT}.*
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/numeric_chrs/chrtruth.vcf \
			 	${DIR}/../data/numeric_chrs/query.vcf \
                -f ${DIR}/../data/numeric_chrs/chrfp.bed \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/numeric_chrs/chrtest.fa \
                -V \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../data/numeric_chrs/expected.summary.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.summary.csv ${DIR}/../data/numeric_chrs/expected.summary.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz | grep -v ^# > ${TMP_OUT}.vcf
diff ${TMP_OUT}.vcf ${DIR}/../data/numeric_chrs/chrexpected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../data/numeric_chrs/chrexpected.vcf"
	exit 1
else
    echo "Chr prefixed truth test successful"
    rm -rf ${TMP_OUT}.*
fi

