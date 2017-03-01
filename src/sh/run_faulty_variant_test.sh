#!/bin/bash

# Test if FP regions are processed accurately
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Test for reading and detecting problematic records: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/open_indel/test.vcf \
			 	${DIR}/../data/open_indel/test_q.vcf \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/open_indel/test.fa -l chrQ \
                -V \
			 	--force-interactive
                # --verbose \

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz | grep -v ^# > ${TMP_OUT}.vcf
diff ${TMP_OUT}.vcf ${DIR}/../data/open_indel/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../data/open_indel/expected.vcf"
	exit 1
else
    echo "Faulty variant test successful"
    rm -rf ${TMP_OUT}.*
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../data/open_indel/test.vcf \
			 	${DIR}/../data/open_indel/test_q_failure.vcf \
			 	-o ${TMP_OUT} \
			 	-X --reference ${DIR}/../data/open_indel/test.fa -l chrQ \
                -V \
			 	--force-interactive
                # --verbose \

if [[ $? != 0 ]]; then
	echo "SUCCESS: faulty variants made hap.py fail."
    rm -rf ${TMP_OUT}.*
else
    echo "FAILURE: hap.py didn't detect faulty variants."
	exit 1
fi

# this trips up pre.py (variant beyond the end of chrQ)
${PYTHON} ${HCDIR}/pre.py \
			 	${DIR}/../data/faulty.vcf \
			 	${TMP_OUT} \
			 	--reference ${DIR}/../data/chrQ.fa -l chrQ \
                --decompose \
			 	--force-interactive
                # --verbose \

if [[ $? != 0 ]]; then
	echo "SUCCESS: faulty variants made pre.py fail."
    rm -rf ${TMP_OUT}.*
else
    echo "FAILURE: pre.py didn't detect faulty variants."
	exit 1
fi

