#!/bin/bash

# Test if FP regions are processed accurately
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Testing FP region matching with default behavior: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/fp_region_accuracy/truth.vcf \
			 	${DIR}/../data/fp_region_accuracy/query.vcf \
                -f ${DIR}/../data/fp_region_accuracy/fp.bed \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/fp_region_accuracy/test.fa -l chrQ \
                -V \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv  ${DIR}/../data/fp_region_accuracy/expected.summary.csv
if [[ $? != 0 ]]; then
	echo "Summary differs! ${TMP_OUT}.summary.csv ${DIR}/../data/fp_region_accuracy/expected.summary.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz | grep -v ^# > ${TMP_OUT}.vcf
diff ${TMP_OUT}.vcf ${DIR}/../data/fp_region_accuracy/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../data/fp_region_accuracy/expected.vcf"
	exit 1
else
    echo "FP region test successful"
    rm -rf ${TMP_OUT}.*
fi

