// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2010-2015 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


/**
 *  \brief Test cases for haplotype set matcher
 *
 *
 * \file test_hapsetmatcher.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <bitset>

#include "HapSetMatcher.hh"
#include "AlleleMatcher.hh"


using namespace variant;

BOOST_AUTO_TEST_CASE(testHapSetMatcher)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                     .parent_path()   // test
                                     .parent_path()   // c++
                                 / boost::filesystem::path("data")
                                 / boost::filesystem::path("chrQ.fa");

    HapSetMatcher hsm(tp.string(), "chrQ");

    hsm.addLeft(RefVar(6, 5, "AAACCC"), 1);
    hsm.addRight(RefVar(0, -1, "AAACCC"), 1);

    HapAssignment assignment;
    hsm.reset(assignment);
    BOOST_CHECK_EQUAL(hsm.optimize(assignment), (unsigned)2l);
    BOOST_CHECK_EQUAL(assignment.variant_assignments.size(), (unsigned)2l);

    auto pass_and_score = hsm.checkAndScore(assignment);
    BOOST_CHECK_EQUAL(pass_and_score.first, true);
    BOOST_CHECK_EQUAL(pass_and_score.second, (unsigned)2l);
    BOOST_CHECK_EQUAL(hsm.numberOfPossibleAssignments(), (unsigned)(3l*3l));
}

BOOST_AUTO_TEST_CASE(testAlleleMatcher)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                     .parent_path()   // test
                                     .parent_path()   // c++
                                 / boost::filesystem::path("data")
                                 / boost::filesystem::path("chrQ.fa");

    AlleleMatcher am(tp.string(), "chrQ");

    am.addLeft(RefVar(6, 5, "AAACCC"), 1);

    // untrimmed allele -- we can trim this to be the same
    am.addRight(RefVar(6, 8, "AAACCCAAA"), 1);

    HapAssignment assignment;
    am.reset(assignment);
    BOOST_CHECK_EQUAL(am.optimize(assignment), (unsigned)2l);
    BOOST_CHECK_EQUAL(assignment.variant_assignments.size(), (unsigned)2l);

    auto pass_and_score = am.checkAndScore(assignment);
    BOOST_CHECK_EQUAL(pass_and_score.first, true);
    BOOST_CHECK_EQUAL(pass_and_score.second, (unsigned)2l);
    BOOST_CHECK_EQUAL(am.numberOfPossibleAssignments(), (unsigned)(2l*2l));
}
