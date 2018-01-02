#!/bin/bash
# Simple consistency test for som.py.
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

#####

echo "# Somatic comparison test for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats.csv!"
	exit 1
fi

#####

echo "# Somatic comparison test with FP regions for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
	exit 1
fi

#####

echo "# Somatic comparison test with FP regions + fixchr_truth for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs_nochr.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
                  --fix-chr-truth \
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
	exit 1
fi

#####

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P --feature-table hcc.strelka.indel"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats.csv!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy_features.py ${TMP_OUT}.features.csv ${DIR}/../../example/sompy/reference_outputs/features.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output features differ:  diff ${TMP_OUT}.features.csv ${DIR}/../../example/sompy/reference_outputs/features.csv"
	exit 1
fi

#####

echo "# Somatic comparison test with FP regions + fixchr_truth + fixchr_query for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_snvs_nochr.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_admix_snvs_nochr.vcf.gz \
                  --fixchr-truth \
                  --fixchr-query \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
if [[ $? != 0 ]]; then
    echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
    exit 1
fi

#####

echo "# Somatic comparison test with FP regions + automatic fixchr_query (no chr prefix in query nor truthset) for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_snvs_nochr.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_admix_snvs_nochr.vcf.gz \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
    exit 1
fi

#####

echo "# Somatic comparison test with FP regions + automatic fixchr_query (chr prefix in query and truthset) for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
    exit 1
fi

#####

echo "# Somatic comparison test with FP regions + automatic fixchr_query (chr prefix in query only) for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_snvs_nochr.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
    exit 1
fi

#####

echo "# Somatic comparison test with FP regions + automatic fixchr_query (chr prefix in truthset only) for ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_admix_snvs_nochr.vcf.gz \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_nochr.bed.gz --count-unk"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_fpr.csv!"
    exit 1
fi

#####

echo "# Somatic comparison test with custom confidence interval leve for ${HCVERSION} from ${HCDIR}l"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_admix_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P --ci-level 0.99"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	  echo "som.py failed!"
	  exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_ci.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	  echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/stats_ci.csv!"
	  exit 1
fi

#####

echo "# Somatic comparison test with FP regions + feature table + happy-stats for ${HCVERSION} from ${HCDIR} - snvs"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_grch38_snvs.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_grch38_admix_pass_snvs.vcf.gz \
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_grch38.bed.gz \
					--feature-table hcc.strelka.snv --count-unk --happy-stats --quiet"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/sompy.strelka_grch38_admix_pass_snvs.stats.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/sompy.strelka_grch38_admix_pass_snvs.stats.csv!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/validate_happy_summary.py --sompy_stats ${TMP_OUT}.stats.csv --happy_summary ${TMP_OUT}.summary.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ between ${TMP_OUT}.stats.csv and ${TMP_OUT}.summary.csv!"
	exit 1
fi

#####

echo "# Somatic comparison test with FP regions + feature table + happy-stats for ${HCVERSION} from ${HCDIR} - indels"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
			 	  ${DIR}/../../example/sompy/PG_admix_truth_grch38_indels.vcf.gz \
			 	  ${DIR}/../../example/sompy/strelka_grch38_admix_pass_indels.vcf.gz \
			 	  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_grch38.bed.gz \
					--feature-table hcc.strelka.indel --count-unk --happy-stats --quiet"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "som.py failed!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/compare_sompy.py ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/sompy.strelka_grch38_admix_pass_indels.stats.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ diff ${TMP_OUT}.stats.csv ${DIR}/../../example/sompy/reference_outputs/sompy.strelka_grch38_admix_pass_indels.stats.csv!"
	exit 1
fi

CMD="${PYTHON} ${DIR}/validate_happy_summary.py --sompy_stats ${TMP_OUT}.stats.csv --happy_summary ${TMP_OUT}.summary.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
	echo "Output counts differ between ${TMP_OUT}.stats.csv and ${TMP_OUT}.summary.csv!"
	exit 1
fi

#####

echo "# Somatic comparison test with FP regions + feature table + AF stratification for ${HCVERSION} from ${HCDIR} - snvs"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_grch38_snvs.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_grch38_admix_pass_snvs.vcf.gz \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_grch38.bed.gz \
                  --feature-table hcc.strelka.snv --count-unk --happy-stats --bin-afs --af-binsize 0.1 --af-truth T_AF --quiet"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/validate_happy_extended.py --sompy_stats ${TMP_OUT}.stats.csv --happy_extended ${TMP_OUT}.extended.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
   echo "Output counts differ between ${TMP_OUT}.stats.csv and ${TMP_OUT}.extended.csv!"
   exit 1
fi

#####

echo "# Somatic comparison test with FP regions + feature table + AF stratification for ${HCVERSION} from ${HCDIR} - indels"

TMP_OUT=`mktemp -t sompy.XXXXXXXXXX`

CMD="${PYTHON} ${HCDIR}/som.py \
                  ${DIR}/../../example/sompy/PG_admix_truth_grch38_indels.vcf.gz \
                  ${DIR}/../../example/sompy/strelka_grch38_admix_pass_indels.vcf.gz \
                  -o ${TMP_OUT} -P -f ${DIR}/../../example/sompy/FP_admix_grch38.bed.gz \
                  --feature-table hcc.strelka.indel --count-unk --happy-stats --bin-afs --af-binsize 0.1 --af-truth T_AF --quiet"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
    echo "som.py failed!"
    exit 1
fi

CMD="${PYTHON} ${DIR}/validate_happy_extended.py --sompy_stats ${TMP_OUT}.stats.csv --happy_extended ${TMP_OUT}.extended.csv"
echo $CMD; $CMD
if [[ $? != 0 ]]; then
   echo "Output counts differ between ${TMP_OUT}.stats.csv and ${TMP_OUT}.extended.csv!"
   exit 1
fi

#####

echo "Som.py test done"
