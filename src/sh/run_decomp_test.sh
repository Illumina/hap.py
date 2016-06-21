#!/bin/bash

# Test variant decomposition
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Variant decomposition test ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/decomp/decomp_test.truth.vcf.gz \
			 	${DIR}/../../example/decomp/decomp_test.query.vcf.gz \
			 	-f ${DIR}/../../example/decomp/decomp_test.conf.bed.gz \
			 	-o ${TMP_OUT} \
			 	-X -V \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

diff ${TMP_OUT}.counts.csv ${DIR}/../../example/decomp/expected.counts.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.counts.csv ${DIR}/../../example/decomp/expected.counts.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz > ${TMP_OUT}.vcf
diff -I ^# ${TMP_OUT}.vcf ${DIR}/../../example/decomp/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../../example/decomp/expected.vcf"
	exit 1
fi

rm -rf ${TMP_OUT}.*

echo "Variant decomposition test was successful."
