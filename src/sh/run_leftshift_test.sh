#!/bin/bash

# Test if left-shifting works or breaks haplotypes
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Testing left-shifts: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/leftshifting_example/truth.vcf \
			 	${DIR}/../data/leftshifting_example/query.vcf \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/leftshifting_example/ref.fa -l chrT \
                --preprocess-truth --leftshift \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv  ${DIR}/../data/leftshifting_example/expected.extended.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! ${TMP_OUT}.extended.csv ${DIR}/../data/leftshifting_example/expected.extended.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz | grep -v ^# > ${TMP_OUT}.vcf
diff ${TMP_OUT}.vcf ${DIR}/../data/leftshifting_example/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../data/leftshifting_example/expected.vcf"
	exit 1
else
    echo "Leftshifting test succeeded."
    rm -rf ${TMP_OUT}.*
fi


