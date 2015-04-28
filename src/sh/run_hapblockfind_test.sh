#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

##############################################################
# Test Hapblockfind
##############################################################

echo "Running hapblockfind test"
TF_r="${DIR}/../data/hapblockfind.bed"
TF_e="${DIR}/../data/expected_hapblockfind.bed"

IF=${DIR}/../data/hapblocks.vcf

cat ${IF} | bgzip > "${IF}.gz"
tabix -p vcf "${IF}.gz"

bin/hapblockfind ${IF}.gz:NA12877 ${IF}.gz:NA12878 \
	-o ${TF_r} -r ${DIR}/../data/microhg19.fa -m 1 --verbose=true -l chr1

diff ${TF_r} ${TF_e}

if [ $? -ne 0 ]; then
	echo "hapblockfind test FAILED. You can inspect ${TF_r} for the failed result."
	exit 1
else
	echo "hapblockfind test SUCCEEDED."
	rm ${TF_r}
fi

if [[ -f "$HG19" ]]; then
	echo "Running Hapblockfind test (2)"

	TF_r=`mktemp -t hapblockfind.XXXXXXXXXX` || exit 1
	TF_e="${DIR}/../../example/performance.blocks.exp.bed"

	bin/hapblockfind ${DIR}/../../example/PG_performance.vcf.gz ${DIR}/../../example/performance.varsonly.vcf.gz \
		-o ${TF_r} -r ${HG19} -l chr1 


	diff ${TF_r} ${TF_e}
	if [ $? -ne 0 ]; then
		echo "hapblockfind test (2) FAILED. You can inspect ${TF_r} for the failed result."
		exit 1
	fi

	${PYTHON} ${DIR}/../python/ovc.py ${TF_r}

	if [ $? -ne 0 ]; then
		echo "hapblockfind test FAILED -- overlaps found. You can inspect ${TF_r} for the failed result."
		exit 1
	fi

	TF_X1=`mktemp -t hapblockfind.XXXXXXXXXX` || exit 1
	TF_X2=`mktemp -t hapblockfind.XXXXXXXXXX` || exit 1
	bin/bcftools view -H ${DIR}/../../example/PG_performance.vcf.gz chr1 | grep -v '0|0' | grep -v '0/0' > ${TF_X1}
	bin/bcftools view -H ${DIR}/../../example/PG_performance.vcf.gz -R ${TF_r} chr1 | grep -v '0|0' | grep -v '0/0' > ${TF_X2}

	diff ${TF_X1} ${TF_X2}	

	if [ $? -ne 0 ]; then
		echo "hapblockfind test FAILED -- doesn't cover all variants. You can inspect ${TF_r} for the failed result."
		exit 1
	else
		rm -f ${TF_X1} ${TF_X2}	
	fi

	TF_X1=`mktemp -t hapblockfind.XXXXXXXXXX` || exit 1
	TF_X2=`mktemp -t hapblockfind.XXXXXXXXXX` || exit 1
	bin/bcftools view -H ${DIR}/../../example/performance.varsonly.vcf.gz chr1 | grep -v '0|0' | grep -v '0/0' > ${TF_X1}
	bin/bcftools view -H ${DIR}/../../example/performance.varsonly.vcf.gz -R ${TF_r} chr1 | grep -v '0|0' | grep -v '0/0' > ${TF_X2}

	diff ${TF_X1} ${TF_X2}	

	if [ $? -ne 0 ]; then
		echo "hapblockfind test FAILED -- doesn't cover all variants. You can inspect ${TF_r} for the failed result."
		exit 1
	else
		rm -f ${TF_X1} ${TF_X2}	
		echo "hapblockfind test (2) SUCCEEDED."
		rm ${TF_r}
	fi

else
	echo "No HG19 file found, Hapblockfind test (2) SKIPPED"
fi
