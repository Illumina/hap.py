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
 *  \brief Test cases for popcount / bit ops
 *
 *
 * \file test_popcount.cpp
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

#include "helpers/Popcount.hh"


BOOST_AUTO_TEST_CASE(testPopcount)
{
    BOOST_CHECK_EQUAL(popcount64(3), 2);
    BOOST_CHECK_EQUAL(popcount64(3 << 10), 2);
    BOOST_CHECK_EQUAL(popcount64(255l << 37), 8);
    BOOST_CHECK_EQUAL(popcount64(-1), 64);
}

BOOST_AUTO_TEST_CASE(testPopcountPermutations)
{
    uint64_t threebits = 0x70;
    typedef std::bitset<64> bin;

//    std::cerr << std::hex << threebits << " / " << bin(threebits) << std::endl;
    threebits = next_permutation(threebits);
    BOOST_CHECK_EQUAL(threebits, (uint64_t)0x83);
//    std::cerr << std::hex << threebits << " / " << bin(threebits) << std::endl;
    threebits = next_permutation(threebits);
    BOOST_CHECK_EQUAL(threebits, (uint64_t)0x85);
//    std::cerr << std::hex << threebits << " / " << bin(threebits) << std::endl;
    threebits = next_permutation(threebits);
    BOOST_CHECK_EQUAL(threebits, (uint64_t)0x86);
//    std::cerr << std::hex << threebits << " / " << bin(threebits) << std::endl;

    threebits = 0x07;
    // 0111 ->
    // 1011 = 8+2+1 = 11
    // 1101 = 8+4+1 = 13
    // 1110 = 8+4+2 = 14
    threebits = next_permutation(threebits);
    BOOST_CHECK_EQUAL(threebits, (uint64_t)0xb);
    threebits = next_permutation(threebits);
    BOOST_CHECK_EQUAL(threebits, (uint64_t)0xd);
    threebits = next_permutation(threebits);
    BOOST_CHECK_EQUAL(threebits, (uint64_t)0xe);
}

BOOST_AUTO_TEST_CASE(testPopcountPermutationsRandom)
{
    typedef std::bitset<64> bin;
    for(int n_test = 0; n_test < 1000; ++n_test)
    {
        uint64_t value = rand64();
        int initial_popcount = popcount64(value);

//        std::cerr << "Test: " << n_test << " : " << std::hex << value << " / " << bin(value) << std::endl;

        uint64_t value_copy = value;
        int popcount_test = 0;
        while(value_copy != 0)
        {
            if(value_copy & 1)
            {
                ++popcount_test;
            }
            value_copy >>= 1;
        }

        BOOST_CHECK_EQUAL(popcount_test, initial_popcount);

        for(int n_perm = 0; n_perm < 1000; ++n_perm)
        {
            value = next_permutation(value);
            BOOST_CHECK_EQUAL(popcount64(value), initial_popcount);
//            std::cerr << "Test: " << n_test << " : " << std::hex << value << " / " << bin(value) << std::endl;
        }
    }
}
