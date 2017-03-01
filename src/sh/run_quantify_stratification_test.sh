#!/bin/bash

# Quantification and GA4GH intermediate file format compliance
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Quantify stratified counting test for ${HCVERSION} from ${HCDIR}"

HG38=${DIR}/../../example/happy/hg38.chr21.fa

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py on NA12878 test file
# note this only works if we don't change the confident regions
# based on the truthset
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_hg38-chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878-GATK3-chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_hg38-chr21.bed.gz \
			 	-r ${HG38} \
			 	-o ${TMP_OUT} \
                --stratification ${DIR}/../../example/happy/stratification.tsv \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# This script checks if the summary precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-stratified.summary.csv
if [[ $? != 0 ]]; then
	echo "All summary differs! -- diff ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-stratified.summary.csv"
	exit 1
fi

# This script checks if the extended precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-stratified.extended.csv
if [[ $? != 0 ]]; then
	echo "All extended differs! -- diff ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-stratified.extended.csv"
	exit 1
fi

rm -rf ${TMP_OUT}.*
