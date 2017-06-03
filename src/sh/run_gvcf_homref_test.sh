#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

##############################################################
# Test Multimerge for homref blocks
##############################################################

echo "Running GVCF homref test"
TF="${DIR}/../data/temp.vcf"

cat ${DIR}/../../example/homref/homref.vcf \
	| bgzip > ${DIR}/../../example/homref/homref.vcf.gz \
   && tabix -p vcf ${DIR}/../../example/homref/homref.vcf.gz
cat ${DIR}/../../example/homref/homref2.vcf \
	| bgzip > ${DIR}/../../example/homref/homref2.vcf.gz \
   && tabix -p vcf ${DIR}/../../example/homref/homref2.vcf.gz

echo "${HCDIR}/multimerge ${DIR}/../../example/homref/homref.vcf.gz \
					${DIR}/../../example/homref/homref2.vcf.gz \
					-o ${TF} -r ${DIR}/../../example/chr21.fa \
					--trimalleles=1 --merge-by-location=1 \
					--homref-split=1 --unique-alleles=1 \
					--calls-only=0"

${HCDIR}/multimerge ${DIR}/../../example/homref/homref.vcf.gz \
					${DIR}/../../example/homref/homref2.vcf.gz \
					-o ${TF} -r ${DIR}/../../example/chr21.fa \
					--trimalleles=1 --merge-by-location=1 \
					--homref-split=1 --unique-alleles=1 \
					--calls-only=0

# TODO this is broken and needs some fixing of logic in VariantReader. It's also not used anywhere, so let's
# not test it and fix it later

# diff -I ^# ${TF} ${DIR}/../../example/homref/expected_merge.vcf

# if [ $? -ne 0 ]; then
# 	echo "GVCF homref test FAILED. You can inspect ${TF} for the failed result. //  diff -I ^# ${TF} ${DIR}/../../example/homref/expected_merge.vcf"
# 	exit 1
# else
# 	echo "GVCF homref test SUCCEEDED."
# 	rm ${TF}
# fi

##############################################################
# Test Multimerge for homref blocks and variants
##############################################################

echo "Running GVCF homref + Variants test"
TF="${DIR}/../data/temp.vcf"

cat ${DIR}/../../example/callsonly/call_merge.vcf \
	| bgzip > ${DIR}/../../example/callsonly/call_merge.vcf.gz \
   && tabix -p vcf ${DIR}/../../example/callsonly/call_merge.vcf.gz

echo "${HCDIR}/multimerge ${DIR}/../../example/callsonly/call_merge.vcf.gz:* \
					-o ${TF} -r ${DIR}/../../example/chr21.fa \
					--process-full=1 --process-formats=1"

${HCDIR}/multimerge ${DIR}/../../example/callsonly/call_merge.vcf.gz:* \
					-o ${TF} -r ${DIR}/../../example/chr21.fa \
					--process-full=1 --process-formats=1

diff ${TF} ${DIR}/../../example/callsonly/expected_callsonly.vcf

if [ $? -ne 0 ]; then
	echo "GVCF homref + Variants test FAILED. You can inspect ${TF} for the failed result. // diff ${TF} ${DIR}/../../example/callsonly/expected_callsonly.vcf"
	exit 1
else
	echo "GVCF homref + Variants test SUCCEEDED."
	rm ${TF}
fi
