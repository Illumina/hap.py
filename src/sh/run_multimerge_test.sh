#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

##############################################################
# Test Multimerge
##############################################################

echo "Running Multimerge test (1)"
TF="${DIR}/../data/temp.vcf"
${HCDIR}/multimerge ${DIR}/../data/merge1.vcf.gz:NA12877 ${DIR}/../data/merge2.vcf.gz:NA12878 -o ${TF} -r ${DIR}/../data/microhg19.fa --trimalleles=1 --merge-by-location=1 --unique-alleles=1

diff ${TF} ${DIR}/../data/expected_merge.vcf

if [ $? -ne 0 ]; then
	echo "Multimerge test FAILED. diff ${TF} ${DIR}/../data/expected_merge.vcf"
	exit 1
else
	echo "Multimerge test SUCCEEDED."
	rm ${TF}
fi

echo "Running Multimerge data import test"
TF="${DIR}/../data/temp.vcf"

cat ${DIR}/../data/import_errors.vcf | bgzip > ${DIR}/../data/import_errors.vcf.gz
tabix -p vcf ${DIR}/../data/import_errors.vcf.gz

${HCDIR}/multimerge ${DIR}/../data/import_errors.vcf.gz:NA12877 -o ${TF} -r ${DIR}/../data/chrQ.fa

diff ${TF} ${DIR}/../data/expected_importtest.vcf

if [ $? -ne 0 ]; then
	echo "Multimerge test FAILED.  diff ${TF} ${DIR}/../data/expected_importtest.vcf"
	exit 1
else
	echo "Multimerge test SUCCEEDED."
	rm ${TF}
fi

if [[ -f "$HG19" ]]; then
	echo "Running Multimerge test (2)"

	HF1=${DIR}/../../example/multimerge/hap_alleles_1.vcf
	HF2=${DIR}/../../example/multimerge/hap_alleles_2.vcf

	cat $HF1 | bgzip > $HF1.gz
	tabix -p vcf $HF1.gz

	cat $HF2 | bgzip > $HF2.gz
	tabix -p vcf $HF2.gz

	TF=`mktemp -t multimerge.XXXXXXXXXX`.vcf
	echo "Output is in $TF"
	${HCDIR}/multimerge $HF1.gz $HF2.gz -r $HG19 -o $TF --leftshift=1 --splitalleles=1

	diff -I^# ${TF} ${DIR}/../../example/multimerge/hap_alleles_leftshifted.vcf

	if [ $? -ne 0 ]; then
		cat $TF
		echo "Multimerge test (2) FAILED. diff -I^# ${TF} ${DIR}/../../example/multimerge/hap_alleles_leftshifted.vcf"
		exit 1
	else
		echo "Multimerge test (2) SUCCEEDED."
	fi

	rm -f TF

	echo "Running Multimerge test (3)"
	# TODO this fails when writing BCF -- fix this

	HF=${DIR}/../../example/multimerge/allele_test.vcf

	cat $HF | bgzip > $HF.gz
	tabix -p vcf $HF.gz

	TF=`mktemp -t multimerge.XXXXXXXXXX`.vcf
	echo "Output is in $TF"
	${HCDIR}/multimerge $HF.gz -r $HG19 -o $TF --leftshift=0

	diff -I^# ${TF} ${DIR}/../../example/multimerge/allele_test.vcf

	if [ $? -ne 0 ]; then
		cat $TF
		echo "Multimerge test (3) FAILED."
		exit 1
	else
		echo "Multimerge test (3) SUCCEEDED."
	fi

	rm -f TF

	echo "Running Multimerge test (4)"

	HF=${DIR}/../../example/multimerge/allele_test.vcf

	cat $HF | bgzip > $HF.gz
	tabix -f -p vcf $HF.gz

	TF=`mktemp -t multimerge.XXXXXXXXXX`.vcf
	echo "Output is in $TF"
	${HCDIR}/multimerge $HF.gz -r $HG19 -o $TF --trimalleles=1 --splitalleles=1 --leftshift=1

	diff -I^# ${TF} ${DIR}/../../example/multimerge/allele_test.sorted.vcf

	if [ $? -ne 0 ]; then
		cat $TF
		echo "Multimerge test (4) FAILED. diff -I^# ${TF} ${DIR}/../../example/multimerge/allele_test.sorted.vcf "
		exit 1
	else
		echo "Multimerge test (4) SUCCEEDED."
	fi

	rm -f TF

	echo "Running Multimerge test (5)"

	HF=${DIR}/../../example/multimerge/features.vcf

	cat $HF | bgzip > $HF.gz
	tabix -f -p vcf $HF.gz

	TF=`mktemp -t multimerge.XXXXXXXXXX`.vcf
	echo "Output is in $TF"
	${HCDIR}/multimerge $HF.gz:* -r $HG19 -o $TF --process-formats=1

	diff -I^# ${TF} ${DIR}/../../example/multimerge/features.processed.vcf

	if [ $? -ne 0 ]; then
		echo "Multimerge test (5) FAILED. diff ${TF} ${DIR}/../../example/multimerge/features.processed.vcf "
		exit 1
	else
		echo "Multimerge test (5) SUCCEEDED."
	fi

	rm -f TF

	# echo "Running Multimerge test (5)"

	# HF=${DIR}/../../example/multimerge/features.vcf

	# cat $HF | bgzip > $HF.gz
	# tabix -f -p vcf $HF.gz

	# TF=`mktemp -t multimerge.XXXXXXXXXX`.vcf
	# echo "Output is in $TF"
	# ${HCDIR}/multimerge $HF.gz:* -r $HG19 -o $TF --process-formats=1 --process-split=1

	# diff -I^# ${TF} ${DIR}/../../example/multimerge/features.processed.split.vcf

	# if [ $? -ne 0 ]; then
	# 	cat $TF
	# 	echo "Multimerge test (6) FAILED. diff -I^# ${TF} ${DIR}/../../example/multimerge/features.processed.split.vcf "
	# 	exit 1
	# else
	# 	echo "Multimerge test (6) SUCCEEDED."
	# fi

	# rm -f TF
else
	echo "Multimerge test (2,3,4, ...) SKIPPED -- need a HG19 fasta file. Please set the environment variable HG19 to point me to that."
fi

exit 0
