// -*- mode: c++; indent-tabs-mode: nil; -*-
//
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
 * \brief Tests for LocationInfo data structure
 *
 *
 * \file test_locationinfo.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "helpers/LocationInfo.hh"

#include <cstdlib>
#include <iostream>

using namespace intervals;

BOOST_AUTO_TEST_CASE(testLocationInfo)
{
    LocationInfo<int> li;

    li.set(100, 10, 20);
    li.set(200, 25, 30);
    li.set(50, 28, 32);

    int s;
    size_t c;
    c = li.querySum(11, 12, s);
    BOOST_CHECK_EQUAL(c, (size_t)2);
    BOOST_CHECK_EQUAL(s, 200);

    c = li.querySum(8, 12, s);
    BOOST_CHECK_EQUAL(c, (size_t)3);
    BOOST_CHECK_EQUAL(s, 300);

    c = li.querySum(18, 22, s);
    BOOST_CHECK_EQUAL(c, (size_t)3);
    BOOST_CHECK_EQUAL(s, 300);

    c = li.querySum(18, 26, s);
    BOOST_CHECK_EQUAL(c, (size_t)5);
    BOOST_CHECK_EQUAL(s, 700);

    c = li.querySum(26, 35, s);
    BOOST_CHECK_EQUAL(c, (size_t)7);
    BOOST_CHECK_EQUAL(s, 2*200 + 5*50);
}

BOOST_AUTO_TEST_CASE(testLocationInfoRandom)
{
    static const int count = 2048;
    static const int icount = 20;
    static const int tcount = 300;

    for (int k = 0; k < tcount; ++k) {
        LocationInfo<int> li;
        int ivs[count];
        bool ivc[count];
        std::fill(std::begin(ivs), std::end(ivs), 0);
        std::fill(std::begin(ivc), std::end(ivc), false);

        for (int i = 0; i < icount; ++i) {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t) count - 1);
            int val = rand() % 10;
            std::fill(ivs + start, ivs + end + 1, val);
            std::fill(ivc + start, ivc + end + 1, true);
            li.set(val, start, end);
        }

        int64_t iv_start = 0;
        int64_t iv_end = count-1;
#ifdef _DEBUG
        // remember for debugging
        LocationInfo<int> li_before = li;
#endif
        if (k % 3 == 0)
        {
            iv_end = (rand()%count);
            if(iv_end > 0)
            {
                li.reset_to(iv_end - 1);
                std::fill(ivc, ivc + iv_end, false);
            }
        }
        else if (k % 3 == 1)
        {
            iv_start = (rand()%count);
            li.reset_from(iv_start);
            std::fill(ivc + iv_start, std::end(ivc), false);
        }

        for (int i = 0; i < count; ++i) {
            int64_t start = rand() % count;
            int64_t orig_start = start;
            int64_t end = std::min(start + rand() % 100, (int64_t) count - 1);

            int sum_li = -1;
            size_t count_li = li.querySum(start, end, sum_li);

            int sum_iv = 0;
            size_t count_iv = 0;
            while(start <= end)
            {
                if (ivc[start])
                {
                    ++count_iv;
                    sum_iv += ivs[start];
                }
                ++start;
            }

            BOOST_CHECK_EQUAL(sum_iv, sum_li);
            BOOST_CHECK_EQUAL(count_iv, count_li);
            if(sum_iv != sum_li || count_iv != count_li)
            {
                std::cerr << "Failed in iteration " << k << " querying from " << orig_start << " to " << end << "\n";
                li.dump(std::cerr);
                return;
            }
        }
    }
}