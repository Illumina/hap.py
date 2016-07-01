#!/bin/bash
#
# Small benchmark using hap.py / vcfeval
#
# Author: Peter Krusche <pkrusche@illumina.com>
#
# This script needs the following things:
#
# * a working hap.py build. This script can e.g. be run from a build folder
# * hap.py must have been built with VCFEval support. This means we need Java
# * There must be a version of Rscript available in the PATH.
# * The ggplot2 package must be installed in R
#
# This script will take a while to run, ca. 20-30min on a 4-core laptop.
# Running with more CPUs will be faster.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
set -e
. ${DIR}/../../src/sh/detect_vars.sh

if [[ -z $HC ]]; then
    HC=$HCDIR/hap.py
fi

if [[ ! -x $HC ]]; then
    echo "Please set HCDIR to point to a hap.py installation."
    exit 1
fi

# These are from the old PG 2.01 release
# PGVCF=$DIR/PG_NA12878_hg38-chr21.vcf.gz
# PGBED=$DIR/PG_Conf_hg38-chr21.bed.gz

# These are from PG-2016.1
PGVCF=$DIR/pg-hg38.bcf
PGBED=$DIR/pg-hg38-conf.bed.gz

REF=$DIR/hg38.chr21.fa

# 1. GATK
# -------
#
# To make ROCs for GATK, we discard the LowQual filter and use QUAL
# For VQSR ROCs, we would use VQSLOD and discard the VQSR Tranche filters

f=GATK3
g=${DIR}/NA12878-GATK3-chr21.vcf.gz
ROC="--roc QUAL --roc-filter LowQual"

echo "----------------------------------------------------------------------------------"
echo
echo "Benchmarking the $f VCF."
# Run four times to compare {xcmp, vcfeval} x {partial credit, no partial credit}
echo
echo "... using xcmp + partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF -o ${f}-xcmp $ROC
echo
echo "... using vcfeval + partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF -o ${f}-ve --engine=vcfeval $ROC
echo
echo "... using xcmp without partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF --no-decompose -o ${f}-nxcmp $ROC
echo
echo "... using vcfeval without partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF --no-decompose -o ${f}-nve --engine=vcfeval $ROC

# 2. Platypus
# -----------
#
# Platypus does not have a simple threshold filter we can remove and draw a ROC curve.
# One way to handle this is to keep all filters and start the ROC at the PASS point,
# which is what we do.

f=Platypus
g=${DIR}/NA12878-Platypus-chr21.vcf.gz
ROC="--roc QUAL"

echo "----------------------------------------------------------------------------------"
echo
echo "Benchmarking the $f VCF."
# Run four times to compare {xcmp, vcfeval} x {partial credit, no partial credit}
echo
echo "... using xcmp + partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF -o ${f}-xcmp $ROC
echo
echo "... using vcfeval + partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF -o ${f}-ve --engine=vcfeval $ROC
echo
echo "... using xcmp without partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF --no-decompose -o ${f}-nxcmp $ROC
echo
echo "... using vcfeval without partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF --no-decompose -o ${f}-nve --engine=vcfeval $ROC

# 3. Freebayes
# ------------
#
# Freebayes has no filters, but we can introduce a filter at QUAL=30 using bcftools

f=Freebayes
g=freebayes-filtered.bcf
ROC="--roc QUAL --roc-filter LowQual"

# pre-filter
$HCDIR/bcftools filter -e 'QUAL<30' -s LowQual ${DIR}/NA12878-Freebayes-chr21.vcf.gz -o freebayes-filtered.bcf -O b

echo "----------------------------------------------------------------------------------"
echo
echo "Benchmarking the $f VCF."
# Run four times to compare {xcmp, vcfeval} x {partial credit, no partial credit}
echo
echo "... using xcmp + partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF -o ${f}-xcmp $ROC
echo
echo "... using vcfeval + partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF -o ${f}-ve --engine=vcfeval $ROC
echo
echo "... using xcmp without partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF --no-decompose -o ${f}-nxcmp $ROC
echo
echo "... using vcfeval without partial credit"
$PYTHON $HC $PGVCF $g -f $PGBED -r $REF --no-decompose -o ${f}-nve --engine=vcfeval $ROC

