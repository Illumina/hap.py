#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

##############################################################
# Test Hapcmp
##############################################################

echo "Running hapcmp test"
TF_r="${DIR}/../data/temp_hapcmp.bed"
TF_e="${DIR}/../data/expected_hapcmp.bed"

if [[ -f "$HG19" ]]; then
	ID="${DIR}/../../example/"
	H_VERSION=$(bin/hapcmp --version)
	${HCDIR}/hapcmp -r $HG19 \
		$ID/hc.bed \
		$ID/hc.vcf.gz \
		$ID/PG_hc.vcf.gz \
		--progress=0 -n 512 \
		--output-sequences=1 \
		--do-alignment=1 \
		-b ${TF_r}

	diff ${TF_r} ${TF_e}

	if [ $? -ne 0 ]; then
		echo "hapcmp test FAILED. See ${TF_r}."
		exit 1
	else
		echo "hapcmp test SUCCEEDED."
		rm ${TF_r}
	fi
else
	echo "hapcmp test SKIPPED. Set the HG19 environment variable to point to a hg19 reference with '>chr21 ...' naming."
fi
