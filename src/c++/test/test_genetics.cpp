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
 *
 * \file test_genetics.cpp
 * \author Richard Shaw
 * \email rshaw@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "helpers/Genetics.hh"

using namespace genetics;

BOOST_AUTO_TEST_CASE(testIsRealBase)
{
    BOOST_CHECK(isRealBase('A') == true);
    BOOST_CHECK(isRealBase('C') == true);
    BOOST_CHECK(isRealBase('G') == true);
    BOOST_CHECK(isRealBase('T') == true);
    BOOST_CHECK(isRealBase('N') == false);
    BOOST_CHECK(isRealBase('D') == false);
    BOOST_CHECK(isRealBase('E') == false);
}

BOOST_AUTO_TEST_CASE(testTransitionBase)
{
    BOOST_CHECK(transitionBase('A') == 'G');
    BOOST_CHECK(transitionBase('C') == 'T');
    BOOST_CHECK(transitionBase('G') == 'A');
    BOOST_CHECK(transitionBase('T') == 'C');
}

BOOST_AUTO_TEST_CASE(testSnvIsTransversion)
{
    bool isValidSnv(false);

    BOOST_CHECK(snvIsTransversion('A', 'C', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('A', 'T', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('A', 'G', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('A', 'N', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == false);


    BOOST_CHECK(snvIsTransversion('C', 'A', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('C', 'G', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('C', 'T', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('C', 'N', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == false);


    BOOST_CHECK(snvIsTransversion('G', 'C', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('G', 'T', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('G', 'A', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('G', 'N', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == false);


    BOOST_CHECK(snvIsTransversion('T', 'A', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('T', 'G', isValidSnv) == true);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('T', 'C', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == true);

    BOOST_CHECK(snvIsTransversion('T', 'N', isValidSnv) == false);
    BOOST_CHECK(isValidSnv == false);
}
