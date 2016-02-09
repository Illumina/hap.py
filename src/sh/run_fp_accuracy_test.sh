#!/bin/bash

# Test if FP regions are processed accurately
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Test for path traversals: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/fp_region_accuracy/truth.vcf \
			 	${DIR}/../data/fp_region_accuracy/query.vcf \
                -f ${DIR}/../data/fp_region_accuracy/fp.bed \
			 	-o ${TMP_OUT} -P \
			 	-X --reference ${DIR}/../data/fp_region_accuracy/test.fa -l chrQ \
                -V \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

diff ${TMP_OUT}.counts.csv ${DIR}/../data/fp_region_accuracy/expected.counts.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.counts.csv ${DIR}/../data/fp_region_accuracy/expected.counts.csv"
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

