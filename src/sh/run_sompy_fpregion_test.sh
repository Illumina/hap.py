#!/bin/bash
# Test for correct handling of FP regions in som.py
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Somatic comparison test for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

# run som.py -- whole "genome" is test region (7BP of chrQ)
${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../data/fp_region_accuracy/truth.vcf \
			 	  ${DIR}/../data/fp_region_accuracy/query.vcf \
                  --reference ${DIR}/../data/fp_region_accuracy/test.fa \
			 	  -o ${TMP_OUT} -P

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

EXPECTED=${DIR}/../data/fp_region_accuracy/expected.sompy.1.csv
${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${EXPECTED}
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${EXPECTED} !"
	exit 1
fi

rm -f ${TMP_OUT}*

# run som.py use FP regions
${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../data/fp_region_accuracy/truth.vcf \
			 	  ${DIR}/../data/fp_region_accuracy/query.vcf \
                  --reference ${DIR}/../data/fp_region_accuracy/test.fa \
                  -f ${DIR}/../data/fp_region_accuracy/fp.bed \
                  --count-unk \
			 	  -o ${TMP_OUT} -P

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

EXPECTED=${DIR}/../data/fp_region_accuracy/expected.sompy.2.csv
${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${EXPECTED}
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${EXPECTED} !"
	exit 1
fi

rm -f ${TMP_OUT}*

# run som.py use FP regions and length
${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../data/fp_region_accuracy/truth.vcf \
			 	  ${DIR}/../data/fp_region_accuracy/query.vcf \
                  --reference ${DIR}/../data/fp_region_accuracy/test.fa \
                  -f ${DIR}/../data/fp_region_accuracy/fp.bed \
                  --count-unk \
                  --fp-region-size 10 \
			 	  -o ${TMP_OUT} -P

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

EXPECTED=${DIR}/../data/fp_region_accuracy/expected.sompy.3.csv
${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${EXPECTED}
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${EXPECTED} !"
	exit 1
fi

rm -f ${TMP_OUT}*

echo "Som.py test done"
