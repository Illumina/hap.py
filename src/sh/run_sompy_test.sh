#!/bin/bash
# Simple consistency test for som.py.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Somatic comparison test for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

# run som.py
${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv !"
	exit 1
fi

echo "Somatic comparison test with FP regions for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

# run som.py
${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix.bed.gz --count-unk

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats_fpr.csv
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats_fpr.csv !"
	exit 1
fi

echo "Somatic comparison test with FP regions + fixchr_truth for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

# run som.py
${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs_nochr.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
                  --fix-chr-truth \
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats_fpr.csv
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats_fpr.csv !"
	exit 1
fi

# run som.py
echo "${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P --feature-table hcc.strelka.indel"

${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P --feature-table hcc.strelka.indel

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv !"
	exit 1
fi

${PYTHON} ${DIR}/compare_sompy_features.py ${TMP_OUT}.features.csv ${DIR}/../../example/sompy/features.csv
if [[ $? != 0 ]]; then
	echo "Output features differ:  diff ${TMP_OUT}.features.csv ${DIR}/../../example/sompy/features.csv"
	exit 1
fi

rm -f ${TMP_OUT}*

echo "Som.py test done"
