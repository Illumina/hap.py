#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

##############################################################
# Test Hapenum
##############################################################

echo "Running hapenum test"
TF="${DIR}/../data/temp.dot"

cat ${DIR}/../data/refgraph1.vcf | bgzip > ${DIR}/../data/refgraph1.vcf.gz && tabix -p vcf -f ${DIR}/../data/refgraph1.vcf.gz
${HCDIR}/hapenum -r ${DIR}/../data/chrQ.fa ${DIR}/../data/refgraph1.vcf.gz:NA12877 \
	--output-dot ${TF} -l chrQ  && dot -Tsvg ${TF} > ${TF}.svg

diff ${TF} ${DIR}/../data/expected_refgraph.dot

if [ $? -ne 0 ]; then
	echo "Hapenum test FAILED. You can inspect ${TF}.svg for the failed result."
	exit 1
else
	echo "Hapenum test SUCCEEDED."
	rm ${TF}
	rm ${TF}.svg
fi
