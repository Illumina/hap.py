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
# Test Hapblockfind
##############################################################

/bin/bash ${DIR}/run_hapblockfind_test.sh

if [[ $? -ne 0 ]]; then
	echo "Hapblockfind test FAILED!"
	exit 1
else
	echo "Hapblockfind test SUCCEEDED!"
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
# Test Hap.py + integration
##############################################################

/bin/bash ${DIR}/run_sompy_test.sh

if [[ $? -ne 0 ]]; then
	echo "Som.py test FAILED!"
	exit 1
else
	echo "Som.py test SUCCEEDED!"
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

