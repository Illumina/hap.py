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

cat ~/workspace_lx/haplotypes/example/homref/homref.vcf \
	| bgzip > ~/workspace_lx/haplotypes/example/homref/homref.vcf.gz \
   && tabix -p vcf ~/workspace_lx/haplotypes/example/homref/homref.vcf.gz
cat ~/workspace_lx/haplotypes/example/homref/homref2.vcf \
	| bgzip > ~/workspace_lx/haplotypes/example/homref/homref2.vcf.gz \
   && tabix -p vcf ~/workspace_lx/haplotypes/example/homref/homref2.vcf.gz

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

diff ${TF} ${DIR}/../../example/homref/expected_merge.vcf

if [ $? -ne 0 ]; then
	echo "GVCF homref test FAILED. You can inspect ${TF} for the failed result."
	exit 1
else
	echo "GVCF homref test SUCCEEDED."
	rm ${TF}
fi

##############################################################
# Test Multimerge for homref blocks and variants
##############################################################

echo "Running GVCF homref + Variants test"
TF="${DIR}/../data/temp.vcf"

cat ~/workspace_lx/haplotypes/example/callsonly/call_merge.vcf \
	| bgzip > ~/workspace_lx/haplotypes/example/callsonly/call_merge.vcf.gz \
   && tabix -p vcf ~/workspace_lx/haplotypes/example/callsonly/call_merge.vcf.gz

echo "${HCDIR}/multimerge ${DIR}/../../example/callsonly/call_merge.vcf.gz:* \
					-o ${TF} -r ${DIR}/../../example/chr21.fa \
					--process-full=1 --process-formats=1"

${HCDIR}/multimerge ${DIR}/../../example/callsonly/call_merge.vcf.gz:* \
					-o ${TF} -r ${DIR}/../../example/chr21.fa \
					--process-full=1 --process-formats=1

diff ${TF} ${DIR}/../../example/callsonly/expected_callsonly.vcf

if [ $? -ne 0 ]; then
	echo "GVCF homref + Variants test FAILED. You can inspect ${TF} for the failed result."
	exit 1
else
	echo "GVCF homref + Variants test SUCCEEDED."
	rm ${TF}
fi
