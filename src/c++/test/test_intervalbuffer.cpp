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
 * \brief Test cases for interval buffers
 *
 *
 * \file test_refvar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "helpers/IntervalBuffer.hh"

#include <cstdlib>
#include <iostream>

using namespace intervals;

BOOST_AUTO_TEST_CASE(testIntervalBuffer)
{
    IntervalBuffer ib;

    ib.addInterval(10, 20, 0);
    ib.addInterval(12, 30, 0);
    ib.addInterval(10, 30, 1);
    ib.addInterval(32, 35, 1);
    ib.addInterval(36, 37, 1);
    ib.addInterval(38, 40, 1);
    ib.addInterval(42, 45, 1);

    IntervalBuffer ib2 = ib;

    BOOST_ASSERT(ib2.isCovered(15, 16, 0));
    BOOST_ASSERT(ib2.isCovered(15, 21, 0));
    BOOST_ASSERT(ib2.isCovered(11, 21, 0));
    BOOST_ASSERT(!ib2.isCovered(11, 31, 0));
    BOOST_ASSERT(!ib2.isCovered(8, 15, 0));
    BOOST_ASSERT(!ib2.isCovered(8, 9, 0));

    BOOST_ASSERT(ib2.isCovered(15, 16, 1));
    BOOST_ASSERT(ib2.isCovered(32, 39, 1));
    BOOST_ASSERT(!ib2.isCovered(32, 43, 1));

    ib.advance(30);

    BOOST_ASSERT(!ib.isCovered(10, 11, 0));
    BOOST_ASSERT(!ib.isCovered(15, 16, 0));
    BOOST_ASSERT(!ib.isCovered(15, 21, 0));
    BOOST_ASSERT(!ib.isCovered(11, 21, 0));
    BOOST_ASSERT(ib.isCovered(30, 30, 0));
    BOOST_ASSERT(!ib.isCovered(8, 15, 0));
    BOOST_ASSERT(!ib.isCovered(8, 9, 0));

    BOOST_ASSERT(!ib.isCovered(15, 16, 1));
    BOOST_ASSERT(ib.isCovered(32, 39, 1));
    BOOST_ASSERT(!ib.isCovered(32, 43, 1));
}

BOOST_AUTO_TEST_CASE(testIntervalBufferRandom)
{
    static const int count = 2048;
    static const int icount = 20;
    static const int tcount = 100;

    for (int k = 0; k < tcount; ++k) {
        bool ivs[count];

        std::vector<std::pair<int64_t, int64_t> > ivlist;

        for (int i = 0; i < count; ++i) {
            ivs[i] = false;
        }

        for (int i = 0; i < icount; ++i) {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t) count - 1);
            ivlist.push_back(std::pair<int64_t, int64_t>(start, end));
            for (int j = start; j <= end; ++j) {
                ivs[j] = true;
            }
        }

        struct iv_less {
            inline bool operator()(const std::pair<int64_t, int64_t> &p,
                                   const std::pair<int64_t, int64_t> &q) {
                return p.first < q.first;
            }
        };

        std::sort(ivlist.begin(), ivlist.end(), iv_less());

        IntervalBuffer ib;

        for (auto const &p : ivlist) {
            ib.addInterval(p.first, p.second, 2);
        }

        for (int i = 0; i < count; ++i) {
            int start = rand() % count;
            int end = std::min(start + rand() % 50, (int) count - 1);

            bool is_covered = true;
            for (int j = start; j <= end; ++j) {
                if (!ivs[j]) {
                    is_covered = false;
                    break;
                }
            }

            if (is_covered != ib.isCovered(start, end, 2)) {
                std::cerr << "(c) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says " <<
                                                                                                     ib.isCovered(start, end, 2) << "\n";
            }

            BOOST_CHECK_EQUAL(is_covered, ib.isCovered(start, end, 2));
        }

        // check overlaps
        for (int i = 0; i < count; ++i) {
            int start = rand() % count;
            int end = std::min(start + rand() % 100, (int) count - 1);

            bool is_covered = false;
            for (int j = start; j <= end; ++j) {
                if (ivs[j]) {
                    is_covered = true;
                    break;
                }
            }

            if (is_covered != ib.hasOverlap(start, end, 2)) {
                std::cerr << "(o) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says " <<
                                                                                                     ib.hasOverlap(start, end, 2) << "\n";
            }

            BOOST_CHECK_EQUAL(is_covered, ib.hasOverlap(start, end, 2));
        }

#ifdef _DEBUG
        // remember for debugging
        IntervalBuffer ib_before = ib;
#endif
        ib.advance(count / 2);

        for (int i = 0; i < count / 2; ++i) {
            ivs[i] = false;
        }

        for (int i = 0; i < count; ++i) {
            int start = rand() % count;
            int end = std::min(start + rand() % 50, (int) count - 1);

            bool is_covered = true;
            for (int j = start; j <= end; ++j) {
                if (!ivs[j]) {
                    is_covered = false;
                    break;
                }
            }

            if (is_covered != ib.isCovered(start, end, 2)) {
                std::cerr << "(c) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says " <<
                                                                                                     ib.isCovered(start, end, 2) << "\n";
            }

            BOOST_CHECK_EQUAL(is_covered, ib.isCovered(start, end, 2));
        }

        // check overlaps
        for (int i = 0; i < count; ++i) {
            int start = rand() % count;
            int end = std::min(start + rand() % 100, (int) count - 1);

            bool is_covered = false;
            for (int j = start; j <= end; ++j) {
                if (ivs[j]) {
                    is_covered = true;
                    break;
                }
            }

            if (is_covered != ib.hasOverlap(start, end, 2)) {
                std::cerr << "(o) Interval " << start << "-" << end << ": ivs say " << is_covered << " ib says " <<
                                                                                                     ib.hasOverlap(start, end, 2) << "\n";
            }

            BOOST_CHECK_EQUAL(is_covered, ib.hasOverlap(start, end, 2));
        }
    }
}

BOOST_AUTO_TEST_CASE(testIntervalBufferOverlap)
{
    IntervalBuffer ib;

    // overlap: 20 - 11 + 1 = 10
    ib.addInterval(10, 20, 0);
    ib.addInterval(11, 30, 1);

    // overlap:   1
    ib.addInterval(110, 110, 0);
    ib.addInterval(100, 110, 1);

    // overlap:   6
    ib.addInterval(120, 125, 0);
    ib.addInterval(120, 125, 1);

    BOOST_CHECK_EQUAL((size_t)17, ib.intersectLanes(0, 1));
}

BOOST_AUTO_TEST_CASE(testIntervalBufferOverlapRandom)
{
    static const int64_t SIZE = 1024;
    srand(100);

    for (int IVL = 1; IVL <= 64; IVL += 3)
    {
        int IVS = SIZE / IVL * 2 / 3;
        for(int trial = 0; trial < 100; ++trial)
        {
            IntervalBuffer ib;
            uint8_t covered[SIZE];

            memset(covered, 0, SIZE);

            for(int q = 0; q < IVS; ++q)
            {
                int mask = 1;
                for(size_t lane = 0; lane < 2; ++lane)
                {
                    const int64_t start = rand() % SIZE;
                    const int64_t end = std::min(SIZE - 1, start + (rand() % IVL) - 1);

                    for(int64_t s = start; s <= end; ++s)
                    {
                        covered[s] |= mask;
                    }

                    ib.addInterval(start, end, lane);
                    mask <<= 1;
                }
            }

            size_t expected = 0;
            for(int64_t p = 0; p < SIZE; ++p)
            {
                if(covered[p] == 3)
                {
                    ++expected;
                }
            }

            const size_t observed1 = ib.intersectLanes(0, 1);
            BOOST_CHECK_EQUAL(expected, observed1);

            const size_t observed2 = ib.intersectLanes(1, 0);
            BOOST_CHECK_EQUAL(expected, observed2);
        }
    }
}
