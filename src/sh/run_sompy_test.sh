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

diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv
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
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix.bed.gz

if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats_fpr.csv
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

diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv
if [[ $? != 0 ]]; then
	echo "Output counts differ diff  ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/stats.csv !"
	exit 1
fi

# HAP-35 : Python < 2.7.8 prints floats differently, so this test won't work
if [[ $(${PYTHON} -c 'import sys; print "0" if tuple(sys.version_info[0:3]) >= (2,7,8) else "1"') == "0" ]];
then
    diff  ${TMP_OUT}.features.csv ${DIR}/../../example/sompy/features.csv
    if [[ $? != 0 ]]; then
    	echo "Output features differ -- diff ${TMP_OUT}.features.csv ${DIR}/../../example/sompy/features.csv !"
    	exit 1
    fi
else
    echo "Python is too old for testing whether the feature tables are the same. Use Python >= 2.7.8"
fi

rm -f ${TMP_OUT}*

echo "Som.py test done"
