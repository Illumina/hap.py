#!/bin/bash

# Simple performance and consistency test.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "PG evaluation test for ${HCVERSION} from ${HCDIR}"

HG19=${DIR}/../../example/chr21.fa

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT} \
			 	-X -P \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# Disabled after v0.2.1-public -- the counts that come out are approximately the same
# but will vary from system to system because allele sorting is not deterministic in this
# version, which leads to slightly different postprocessing.
#
# diff -I fileDate -I source_version ${TMP_OUT}.counts.json ${DIR}/../../example/happy/expected.counts.json
# if [[ $? != 0 ]]; then
# 	echo "Counts differ! -- diff ${TMP_OUT}.counts.json ${DIR}/../../example/happy/expected.counts.json"
# 	exit 1
# fi

# This script checks if the summary precision / recall figures have changed significantly
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/happy/expected.summary.csv
if [[ $? != 0 ]]; then
	echo "PASS summary differs! -- diff ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/happy/expected.pass.summary.csv"
	exit 1
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT}.pass \
			 	-X \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# Disabled after v0.2.1-public -- the counts that come out are approximately the same
# but will vary from system to system because allele sorting is not deterministic in this
# version, which leads to slightly different postprocessing.
#
# diff -I fileDate -I source_version ${TMP_OUT}.pass.counts.json ${DIR}/../../example/happy/expected.pass.counts.json
# if [[ $? != 0 ]]; then
# 	echo "PASS counts differ! -- diff ${TMP_OUT}.pass.counts.json ${DIR}/../../example/happy/expected.pass.counts.json"
# 	exit 1
# fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/happy/expected.pass.summary.csv
if [[ $? != 0 ]]; then
	echo "PASS summary differs! -- diff ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/happy/expected.pass.summary.csv"
	exit 1
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/happy/PG_NA12878_chr21.vcf.gz \
			 	${DIR}/../../example/happy/NA12878_chr21.vcf.gz \
			 	-f ${DIR}/../../example/happy/PG_Conf_chr21.bed.gz \
			 	-r ${DIR}/../../example/chr21.fa \
			 	-o ${TMP_OUT}.unhappy \
			 	-X --unhappy \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# Disabled after v0.2.1-public -- the counts that come out are approximately the same
# but will vary from system to system because allele sorting is not deterministic in this
# version, which leads to slightly different postprocessing.
#
# diff -I fileDate -I source_version ${TMP_OUT}.pass.counts.json ${DIR}/../../example/happy/expected.pass.counts.json
# if [[ $? != 0 ]]; then
# 	echo "PASS counts differ! -- diff ${TMP_OUT}.pass.counts.json ${DIR}/../../example/happy/expected.pass.counts.json"
# 	exit 1
# fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.unhappy.summary.csv ${DIR}/../../example/happy/expected.unhappy.summary.csv
if [[ $? != 0 ]]; then
	echo "PASS summary differs! -- diff ${TMP_OUT}.unhappy.summary.csv ${DIR}/../../example/happy/expected.unhappy.summary.csv"
	exit 1
fi

rm -f ${TMP_OUT}*
echo "PG quantification test SUCCEEDED!"

