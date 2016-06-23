#!/bin/bash

# Simple performance and consistency test.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Integration test for ${HCVERSION} from ${HCDIR}"

cat ${DIR}/../../example/integration/integrationtest_lhs.vcf | bgzip > ${DIR}/../../example/integration/integrationtest_lhs.vcf.gz
cat ${DIR}/../../example/integration/integrationtest_rhs.vcf | bgzip > ${DIR}/../../example/integration/integrationtest_rhs.vcf.gz

tabix -f -p vcf ${DIR}/../../example/integration/integrationtest_lhs.vcf.gz
tabix -f -p vcf ${DIR}/../../example/integration/integrationtest_rhs.vcf.gz

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# fallback HG19 locations
if [[ ! -f "$HG19" ]]; then
	HG19=/illumina/development/iSAAC/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
fi

if [[ ! -f "$HG19" ]]; then
	HG19=~/workspace/human_genome/hg19.fa
fi

if [[ -f "$HG19" ]]; then
	TMP_OUT0=`mktemp -t multimerge.XXXXXXXXXX`
	# multimerge and compare
	${HCDIR}/multimerge ${DIR}/../../example/integration/integrationtest_lhs.vcf.gz  ${DIR}/../../example/integration/integrationtest_rhs.vcf.gz \
		-o $TMP_OUT0.vcf -r $HG19 --process-full 1

	diff -I '^#' ${TMP_OUT0}.vcf ${DIR}/../../example/integration/integrationtest_merged.vcf
	if [[ $? != 0 ]]; then
		echo "Merged integration VCF doesn't match -- inspect ${TMP_OUT0}.vcf vs ${DIR}/../../example/integration/integrationtest_merged.vcf."
		exit 1
	else
		rm -f ${TMP_OUT0}.vcf
	fi
else
	echo "Multimerge integration test SKIPPED (set up the HG19 environment variable to point to the reference fasta file)"
	exit 1
fi

# HAP-143 run hap.py with one empty input
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/integration/integrationtest_empty.vcf \
			 	${DIR}/../../example/integration/integrationtest_rhs.vcf.gz \
			 	-o ${TMP_OUT}.e0 \
			 	-V -X --output-vtc \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/integration/integrationtest_lhs.vcf.gz \
			 	${DIR}/../../example/integration/integrationtest_empty.vcf \
			 	-o ${TMP_OUT}.e1 \
			 	-V -X --output-vtc \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/integration/integrationtest_lhs.vcf.gz \
			 	${DIR}/../../example/integration/integrationtest_rhs.vcf.gz \
			 	-o ${TMP_OUT} \
			 	-V -X --output-vtc \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

cat ${TMP_OUT}.vcf.gz | gunzip | grep -v ^# > ${TMP_OUT}.vcf

# run hap.py, without any clever comparison
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/integration/integrationtest_lhs.vcf.gz \
			 	${DIR}/../../example/integration/integrationtest_rhs.vcf.gz \
			 	-o ${TMP_OUT}.unhappy \
			 	-V -X --output-vtc \
			 	--force-interactive --unhappy

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

cat ${TMP_OUT}.unhappy.vcf.gz | gunzip | grep -v ^# > ${TMP_OUT}.unhappy.vcf

# run hap.py for PASS only, without VTC
${PYTHON} ${HCDIR}/hap.py \
			 	-l chr21 \
			 	${DIR}/../../example/integration/integrationtest_lhs.vcf.gz \
			 	${DIR}/../../example/integration/integrationtest_rhs.vcf.gz \
			 	-o ${TMP_OUT}.pass \
			 	-V -X --pass-only \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

cat ${TMP_OUT}.pass.vcf.gz | gunzip | grep -v ^# > ${TMP_OUT}.pass.vcf

diff -I fileDate -I source_version ${TMP_OUT}.vcf ${DIR}/../../example/integration/integrationtest.vcf
if [[ $? != 0 ]]; then
	echo "Output variants differ -- vimdiff ${TMP_OUT}.vcf ${DIR}/../../example/integration/integrationtest.vcf !"
	exit 1
fi

diff -I fileDate -I source_version ${TMP_OUT}.unhappy.vcf ${DIR}/../../example/integration/integrationtest.unhappy.vcf
if [[ $? != 0 ]]; then
	echo "Output variants differ -- vimdiff ${TMP_OUT}.unhappy.vcf ${DIR}/../../example/integration/integrationtest.unhappy.vcf !"
	exit 1
fi

diff -I fileDate -I source_version ${TMP_OUT}.pass.vcf ${DIR}/../../example/integration/integrationtest.pass.vcf
if [[ $? != 0 ]]; then
	echo "Pass output variants differ -- vimdiff ${TMP_OUT}.pass.vcf ${DIR}/../../example/integration/integrationtest.pass.vcf !"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/integration/integrationtest.summary.csv
if [[ $? != 0 ]]; then
	echo "Summary differs! ${TMP_OUT}.summary.csv ${DIR}/../../example/integration/integrationtest.summary.csv"
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/integration/integrationtest.summary.pass.csv
if [[ $? != 0 ]]; then
	echo "Pass summary differs! ${TMP_OUT}.pass.summary.csv ${DIR}/../../example/integration/integrationtest.summary.pass.csv"
	exit 1
fi

## performance test. this does half of chr21

# run hap.py
${PYTHON} ${HCDIR}/hap.py $@ \
			 	-l chr21 \
			 	${DIR}/../../example/PG_performance.vcf.gz \
			 	${DIR}/../../example/performance.vcf.gz \
			 	-o ${TMP_OUT}.performance \
			 	-V -X --output-vtc \
			 	-f ${DIR}/../../example/performance.confident.bed.gz \
                -r ${DIR}/../../example/chr21.fa \
			 	--force-interactive --threads 4

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

cat ${TMP_OUT}.performance.vcf.gz | gunzip | grep -v ^# > ${TMP_OUT}.performance.vcf

diff -I ^# ${TMP_OUT}.performance.vcf ${DIR}/../../example/integration/integrationtest.performance.vcf
if [[ $? != 0 ]]; then
	echo "Performance output variants differ! vimdiff ${TMP_OUT}.performance.vcf ${DIR}/../../example/integration/integrationtest.performance.vcf "
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.performance.summary.csv ${DIR}/../../example/integration/integrationtest.performance.summary.csv
if [[ $? != 0 ]]; then
	echo "Pass summary differs! vimdiff ${TMP_OUT}.performance.summary.csv ${DIR}/../../example/integration/integrationtest.performance.summary.csv"
	exit 1
fi

# single-threaded version should produce same result
# run hap.py
${PYTHON} ${HCDIR}/hap.py $@ \
			 	-l chr21 \
			 	${DIR}/../../example/PG_performance.vcf.gz \
			 	${DIR}/../../example/performance.vcf.gz \
			 	-o ${TMP_OUT}.performance.t1 \
			 	-V -X --output-vtc \
			 	-f ${DIR}/../../example/performance.confident.bed.gz \
                -r ${DIR}/../../example/chr21.fa \
			 	--force-interactive --threads 1

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

cat ${TMP_OUT}.performance.t1.vcf.gz | gunzip | grep -v ^# > ${TMP_OUT}.performance.t1.vcf

echo "Comparing vcfs"
diff -I ^# ${TMP_OUT}.performance.t1.vcf ${DIR}/../../example/integration/integrationtest.performance.vcf
if [[ $? != 0 ]]; then
    echo "Performance output variants differ (t1)! vimdiff ${TMP_OUT}.performance.t1.vcf ${DIR}/../../example/integration/integrationtest.performance.vcf "
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.performance.t1.summary.csv ${DIR}/../../example/integration/integrationtest.performance.summary.csv
if [[ $? != 0 ]]; then
	echo "Pass summary differs! vimdiff ${TMP_OUT}.performance.t1.summary.csv ${DIR}/../../example/integration/integrationtest.performance.summary.csv"
	exit 1
fi

# BCF version should produce same result
# run hap.py
${PYTHON} ${HCDIR}/hap.py $@ \
			 	-l chr21 \
			 	${DIR}/../../example/PG_performance.vcf.gz \
			 	${DIR}/../../example/performance.vcf.gz \
			 	-o ${TMP_OUT}.performance.bcf \
			 	-V -X --output-vtc --bcf \
			 	-f ${DIR}/../../example/performance.confident.bed.gz \
                -r ${DIR}/../../example/chr21.fa \
			 	--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

${HCDIR}/bcftools view -H ${TMP_OUT}.performance.bcf.bcf -o ${TMP_OUT}.performance.bcf.vcf -O v

echo "Comparing vcfs"
diff -I ^# ${TMP_OUT}.performance.bcf.vcf ${DIR}/../../example/integration/integrationtest.performance.vcf
if [[ $? != 0 ]]; then
    echo "Performance output variants differ (bcf)! vimdiff ${TMP_OUT}.performance.bcf.vcf ${DIR}/../../example/integration/integrationtest.performance.vcf "
	exit 1
fi

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.performance.bcf.summary.csv ${DIR}/../../example/integration/integrationtest.performance.summary.csv
if [[ $? != 0 ]]; then
	echo "Pass summary differs! vimdiff ${TMP_OUT}.performance.bcf.summary.csv ${DIR}/../../example/integration/integrationtest.performance.summary.csv"
	exit 1
fi

rm -f ${TMP_OUT}*
echo "Integration test SUCCEEDED!"

