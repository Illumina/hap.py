#!/bin/bash

##############################################################
# Test setup
##############################################################

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. ${DIR}/detect_vars.sh

##############################################################
# Test blocksplit
##############################################################


echo "Running blocksplit test"

TF_r=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

VCF1=${DIR}/../../example/happy/PG_NA12878_hg38-chr21.vcf.gz
VCF2=${DIR}/../../example/happy/NA12878-GATK3-chr21.vcf.gz

${HCDIR}/blocksplit $VCF1 $VCF2 \
	-o ${TF_r} -l chr21 -w 10000

${PYTHON} ${DIR}/../python/ovc.py ${TF_r}
if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- overlaps found. You can inspect ${TF_r} for the failed result."
	exit 1
fi

TF_X1=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
TF_X2=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

${HCDIR}/bcftools view -H $VCF1 chr21 > ${TF_X1}
${HCDIR}/bcftools view -H $VCF1 chr21 -R ${TF_r} > ${TF_X2}

diff ${TF_X1} ${TF_X2}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- doesn't cover all variants. You can inspect ${TF_r} / ${TF_X1} / ${TF_X2} for the failed result."
	exit 1
else
	rm -f ${TF_X1} ${TF_X2}
fi

TF_X1=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
TF_X2=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

${HCDIR}/bcftools view -H $VCF2 chr21 > ${TF_X1}
${HCDIR}/bcftools view -H $VCF2 chr21 -R ${TF_r} > ${TF_X2}

diff ${TF_X1} ${TF_X2}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- doesn't cover all variants. You can inspect ${TF_r} / ${TF_X1} / ${TF_X2} for the failed result."
	exit 1
else
	rm -f ${TF_X1} ${TF_X2}
fi

echo "Blocksplit test was successful"
