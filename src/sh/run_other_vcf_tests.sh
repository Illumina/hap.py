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
			 	${DIR}/../data/per_sample_ft_lhs.vcf \
			 	${DIR}/../data/per_sample_ft_rhs.vcf \
			 	-o ${TMP_OUT} \
			 	--reference ${DIR}/../data/chrQ.fa \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../data/per_sample_ft_summary.csv
if [[ $? != 0 ]]; then
	echo "Counts differ! diff ${TMP_OUT}.summary.csv ${DIR}/../data/per_sample_ft_summary.csv"
	exit 1
else
    echo "Variant filtering test successful"
    rm -rf ${TMP_OUT}.*
fi

####################################
# haploid variants (chrX/Y handling)
####################################

echo "Test for handling haploid records: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/haploid/truth.vcf \
			 	${DIR}/../../example/haploid/query.vcf \
			 	-o ${TMP_OUT} \
			 	-r ${DIR}/../../example/haploid/test.fa \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

cat ${TMP_OUT}.vcf.gz | gunzip | grep -v ^# > ${TMP_OUT}.vcf
diff -I '^#' ${TMP_OUT}.vcf ${DIR}/../../example/haploid/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff -I '^#' ${TMP_OUT}.vcf ${DIR}/../../example/haploid/expected.vcf"
	exit 1
else
    echo "Haploid test successful"
    rm -rf ${TMP_OUT}.*
fi


