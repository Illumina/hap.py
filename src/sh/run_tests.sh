#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

##############################################################
# Boost unit tests
##############################################################

# echo "Running BOOST_TEST tests."

bin/test_haplotypes
if [[ $? -ne 0 ]]; then
	echo "Boost Unit tests FAILED!"
	exit 1
else
	echo "Boost Unit tests SUCCEEDED!"
fi

##############################################################
# Test Multimerge
##############################################################

/bin/bash ${DIR}/run_multimerge_test.sh

if [[ $? -ne 0 ]]; then
	echo "Multimerge tests FAILED!"
	exit 1
else
	echo "Multimerge tests SUCCEEDED!"
fi

##############################################################
# Test Hapenum
##############################################################

/bin/bash ${DIR}/run_hapenum_test.sh

if [[ $? -ne 0 ]]; then
	echo "Hapenum test FAILED!"
	exit 1
else
	echo "Hapenum test SUCCEEDED!"
fi

##############################################################
# Test Hapcmp
##############################################################

/bin/bash ${DIR}/run_hapcmp_test.sh

if [[ $? -ne 0 ]]; then
	echo "Hapcmp test FAILED!"
	exit 1
else
	echo "Hapcmp test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + path traversals
##############################################################

/bin/bash ${DIR}/run_pathtraversal_test.sh

if [[ $? -ne 0 ]]; then
	echo "Path traversal test FAILED!"
	exit 1
else
	echo "Path traversal test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + FP regions
##############################################################

/bin/bash ${DIR}/run_fp_accuracy_test.sh

if [[ $? -ne 0 ]]; then
	echo "FP region test FAILED!"
	exit 1
else
	echo "FP region test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + Faulty input variants
##############################################################

/bin/bash ${DIR}/run_faulty_variant_test.sh

if [[ $? -ne 0 ]]; then
	echo "Faulty variant test FAILED!"
	exit 1
else
	echo "Faulty variant test SUCCEEDED!"
fi

##############################################################
# Test Hap.py safe leftshifting
##############################################################

/bin/bash ${DIR}/run_leftshift_test.sh

if [[ $? -ne 0 ]]; then
	echo "Leftshift test FAILED!"
	exit 1
else
	echo "Leftshift test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + other VCF items
##############################################################

/bin/bash ${DIR}/run_other_vcf_tests.sh

if [[ $? -ne 0 ]]; then
	echo "Other VCF tests FAILED!"
	exit 1
else
	echo "Other VCF tests SUCCEEDED!"
fi

##############################################################
# Test hom-ref block expansion and calls-only preprocessing
##############################################################

/bin/bash ${DIR}/run_gvcf_homref_test.sh

if [[ $? -ne 0 ]]; then
	echo "GVCF hom-ref test FAILED!"
	exit 1
else
	echo "GVCF hom-ref test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + chr prefix detection
##############################################################

/bin/bash ${DIR}/run_chrprefix_test.sh

if [[ $? -ne 0 ]]; then
	echo "Chr prefix detection tests FAILED!"
	exit 1
else
	echo "Chr prefix detection tests SUCCEEDED!"
fi

##############################################################
# Test Hap.py variant decomposition into primitives
##############################################################

/bin/bash ${DIR}/run_decomp_test.sh

if [[ $? -ne 0 ]]; then
	echo "Variant decomposition test FAILED!"
	exit 1
else
	echo "Variant decomposition test SUCCEEDED!"
fi

##############################################################
# Test contig length calculation
##############################################################

${PYTHON} ${DIR}/run_fastasize_test.py
if [[ $? -ne 0 ]]; then
    echo "Contig length calculation test FAILED!"
    exit 1
else
    echo "Contig length calculation test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + integration
##############################################################

/bin/bash ${DIR}/run_integration_test.sh

if [[ $? -ne 0 ]]; then
	echo "Integration test FAILED!"
	exit 1
else
	echo "Integration test SUCCEEDED!"
fi

##############################################################
# Test Hap.py scmp allele and distance-based comparison
##############################################################

/bin/bash ${DIR}/run_scmp_test.sh

if [[ $? -ne 0 ]]; then
	echo "Integration test FAILED!"
	exit 1
else
	echo "Integration test SUCCEEDED!"
fi

##############################################################
# Test Hap.py on tricky test cases
##############################################################

/bin/bash ${DIR}/run_giab_test.sh

if [[ $? -ne 0 ]]; then
	echo "Tricky indel test FAILED!"
	exit 1
else
	echo "Tricky indel test SUCCEEDED!"
fi

##############################################################
# Test Performance + Consistency
##############################################################

/bin/bash ${DIR}/run_performance_test.sh

if [[ $? -ne 0 ]]; then
	echo "Performance / Consistency test FAILED!"
	exit 1
else
	echo "Performance / Consistency test SUCCEEDED!"
fi

##############################################################
# Test GA4GH quantification
##############################################################

/bin/bash ${DIR}/run_quantify_test.sh

if [[ $? -ne 0 ]]; then
	echo "Quantify integration test FAILED!"
	exit 1
else
	echo "Quantify integration test SUCCEEDED!"
fi

##############################################################
# Test GA4GH stratified quantification
##############################################################

/bin/bash ${DIR}/run_quantify_stratification_test.sh

if [[ $? -ne 0 ]]; then
	echo "Quantify stratification test FAILED!"
	exit 1
else
	echo "Quantify stratification test SUCCEEDED!"
fi


##############################################################
# Test PG Counting
##############################################################

/bin/bash ${DIR}/run_happy_pg_test.sh

if [[ $? -ne 0 ]]; then
	echo "PG integration test FAILED!"
	exit 1
else
	echo "PG integration test SUCCEEDED!"
fi

##############################################################
# Test Hap.py + integration
##############################################################

/bin/bash ${DIR}/run_sompy_test.sh

if [[ $? -ne 0 ]]; then
	echo "Som.py test FAILED!"
	exit 1
else
	echo "Som.py test SUCCEEDED!"
fi
