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

${HCDIR}/blocksplit ${DIR}/../../example/PG_performance.vcf.gz ${DIR}/../../example/performance.varsonly.vcf.gz \
	-o ${TF_r} -l chr1

${PYTHON} ${DIR}/../python/ovc.py ${TF_r}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- overlaps found. You can inspect ${TF_r} for the failed result."
	exit 1
fi

TF_X1=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
TF_X2=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

${HCDIR}/bcftools view -H ${DIR}/../../example/PG_performance.vcf.gz chr1 > ${TF_X1}
${HCDIR}/bcftools view -H ${DIR}/../../example/PG_performance.vcf.gz -R ${TF_r} > ${TF_X2}

diff ${TF_X1} ${TF_X2}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- doesn't cover all variants. You can inspect ${TF_r} / ${TF_X1} / ${TF_X2} for the failed result."
	exit 1
else
	rm -f ${TF_X1} ${TF_X2}
fi

TF_X1=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
TF_X2=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
${HCDIR}/bcftools view -H ${DIR}/../../example/performance.varsonly.vcf.gz chr1 > ${TF_X1}
${HCDIR}/bcftools view -H ${DIR}/../../example/performance.varsonly.vcf.gz -R ${TF_r} > ${TF_X2}

diff ${TF_X1} ${TF_X2}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- doesn't cover all variants. You can inspect ${TF_r} / ${TF_X1} / ${TF_X2} for the failed result."
	exit 1
else
	rm -f ${TF_X1} ${TF_X2}
fi
