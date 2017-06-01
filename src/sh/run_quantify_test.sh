#!/bin/bash

# Quantification and GA4GH intermediate file format compliance
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Quantify evaluation test for ${HCVERSION} from ${HCDIR}"

HG19=${DIR}/../../example/chr21.fa

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py on NA12878 test file
# note this only works if we don't change the confident regions
# based on the truthset
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT} \
                --no-adjust-conf-regions \
			 	--force-interactive --verbose

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# This script checks if the summary precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-qfy.summary.csv
if [[ $? != 0 ]]; then
	echo "All summary differs! -- diff ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected-qfy.summary.csv"
	exit 1
fi

# This script checks if the extended precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_extended.py ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-qfy.extended.csv
if [[ $? != 0 ]]; then
	echo "All extended differs! -- diff ${TMP_OUT}.extended.csv ${DIR}/../../example/happy/expected-qfy.extended.csv"
	exit 1
fi


# hap.py writes GA4GH-compliant output VCFs
# we should be able to re-quantify using the GA4GH
# spec and get the same result
# Note that the default is now to run with --adjust-conf-regions in hap.py,
# so we have to also pass the truth VCF here
${PYTHON} ${HCDIR}/qfy.py \
           ${TMP_OUT}.vcf.gz \
           -o ${TMP_OUT}.qfy  \
           -r ${DIR}/../../example/chr21.fa \
           -f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
           -t ga4gh -X --verbose

if [[ $? != 0 ]]; then
	echo "qfy.py failed!"
	exit 1
fi

# these are all the metrics. We filter out all the command line bits
gunzip -c ${TMP_OUT}.metrics.json.gz | ${PYTHON} -mjson.tool | grep -v timestamp | grep -v hap.py > ${TMP_OUT}.hap.m.json
if [[ $? != 0 ]] || [[ ! -s ${TMP_OUT}.hap.m.json ]]; then
	echo "Cannot unzip metrics for original run."
	exit 1
fi
gunzip -c ${TMP_OUT}.qfy.metrics.json.gz | ${PYTHON} -mjson.tool | grep -v timestamp | grep -v qfy.py > ${TMP_OUT}.qfy.m.json
if [[ $? != 0 ]] || [[ ! -s ${TMP_OUT}.qfy.m.json ]]; then
	echo "Cannot unzip metrics for re-quantified run."
	exit 1
fi

# This script checks if the summary precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.qfy.summary.csv ${DIR}/../../example/happy/expected-qfy.summary.csv
if [[ $? != 0 ]]; then
	echo "All summary differs! -- diff ${TMP_OUT}.qfy.summary.csv ${DIR}/../../example/happy/expected-qfy.summary.csv"
	exit 1
fi

diff ${TMP_OUT}.hap.m.json ${TMP_OUT}.qfy.m.json > /dev/null
if [[ $? != 0 ]]; then
	echo "Re-quantified counts are different! diff ${TMP_OUT}.hap.m.json ${TMP_OUT}.qfy.m.json "
	exit 1
fi

rm -rf ${TMP_OUT}.*
