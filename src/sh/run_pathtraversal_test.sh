#!/bin/bash

# Test invalid path traversals
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Test for path traversals: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
				${DIR}/../data/pathtraversal/test.vcf \
				${DIR}/../data/pathtraversal/test2.vcf \
				-o ${TMP_OUT} -P \
				-X --reference ${DIR}/../data/pathtraversal/test.fa -l chrQ \
				--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

diff ${TMP_OUT}.counts.csv ${DIR}/../data/pathtraversal/expected.counts.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.counts.csv ${DIR}/../data/pathtraversal/expected.counts.csv"
	exit 1
else
    echo "Path traversal test successful"
    rm -rf ${TMP_OUT}.*
fi

