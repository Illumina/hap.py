#!/bin/bash

# Simple performance and consistency test.
# We compare the output of simplecmp to hapcmp
#
# If all goes well, we don't get any blocks in which
# simplecmp produces a match, but hapcmp produces a mismatch.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

if [[ ! -f "${HCDIR}/xcmp" ]]; then
	echo "Cannot find HC binaries. Set HCDIR to the bin directory."
	exit 1
fi

${HCDIR}/xcmp \
	${DIR}/../../example/PG_performance.vcf.gz \
	${DIR}/../../example/performance.vcf.gz \
	--always-hapcmp 1 -e performance.bed -o performance.vcf \
	-r ${HG19} -f 0 -n 256

if [[ $? -ne 0 ]]; then
	echo "xcmp failed!"
	exit 1
fi

# find differences
SIMPLE_MIS_HC_MATCH=`cat performance.bed | grep 'simple:match' | wc -l`
SIMPLE_MIS_HC_MIS=`cat performance.bed | grep 'mismatch' | wc -l`
SIMPLE_MATCH_HC_MIS=`cat performance.bed | grep 'suspicious' | wc -l`

echo "Mismatches in both: '$SIMPLE_MIS_HC_MIS' (== cases not rescued through haplotype comparison)"
echo "Matches in hapcmp but not simplecmp: '$SIMPLE_MIS_HC_MATCH' (== cases rescued through haplotype comparison)"
echo "Suspicious matches: '$SIMPLE_MATCH_HC_MIS' (if this is >0, things are not good)"

if [[ $SIMPLE_MATCH_HC_MIS > 0 ]]; then
	echo "ERROR: Consistency check failed: there are blocks for which simplecmp reports a match, but hapcmp reports a mismatch (on VCF)."
	exit 1
else
	echo "SUCCESS: Performance + consistency check (VCF)."
fi

## GVCF test

${HCDIR}/xcmp \
	${DIR}/../../example/PG_performance.vcf.gz \
	${DIR}/../../example/performance.vcf.gz \
	-r ${HG19} -f 0 -n 256  \
	--always-hapcmp 1 -e performance.bed -o performance.vcf

if [[ $? -ne 0 ]]; then
	echo "xcmp failed!"
	exit 1
fi

# find differences
SIMPLE_MIS_HC_MATCH=`cat performance.bed | grep 'simple:match' | wc -l`
SIMPLE_MIS_HC_MIS=`cat performance.bed | grep 'mismatch' | wc -l`
SIMPLE_MATCH_HC_MIS=`cat performance.bed | grep 'suspicious' | wc -l`

echo "Mismatches in both: '$SIMPLE_MIS_HC_MIS' (== cases not rescued through haplotype comparison)"
echo "Matches in hapcmp but not simplecmp: '$SIMPLE_MIS_HC_MATCH' (== cases rescued through haplotype comparison)"
echo "Suspicious matches: '$SIMPLE_MATCH_HC_MIS' (if this is >0, things are not good)"

if [[ $SIMPLE_MATCH_HC_MIS > 0 ]]; then
	echo "ERROR: Consistency check failed: there are blocks for which simplecmp reports a match, but hapcmp reports a mismatch (on VCF)."
	exit 1
else
	echo "SUCCESS: Performance + consistency check (VCF)."
fi
