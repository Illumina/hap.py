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

#include "helpers/IntervalList.hh"

#include <cstdlib>
#include <iostream>

using namespace intervals;


template<class _t>
void listcmp(
    std::list<_t> expected,
    std::list<_t> actual)
{
    auto i1 = expected.begin();
    auto i2 = actual.begin();
    size_t count = 0;
    while(i1 != expected.end() && i2 != actual.end())
    {
        BOOST_CHECK_EQUAL(*i1, *i2);
        ++count;
        ++i1;
        ++i2;
    }
    BOOST_CHECK_EQUAL(count, expected.size());
}

template<>
void listcmp(
    std::list<std::pair<int64_t, int64_t> > expected,
    std::list<std::pair<int64_t, int64_t> > actual)
{
    auto i1 = expected.begin();
    auto i2 = actual.begin();
    size_t count = 0;
    while(i1 != expected.end() && i2 != actual.end())
    {
        BOOST_CHECK_EQUAL(i1->first, i2->first);
        BOOST_CHECK_EQUAL(i1->second, i2->second);
        ++count;
        ++i1;
        ++i2;
    }
    BOOST_CHECK_EQUAL(count, expected.size());
}

struct ivcount : public interval
{
    ivcount() : interval(), count(0) {}
    ivcount(int64_t _start, int64_t _end, int _count) : interval(_start, _end), count(_count) {}
    virtual void merge(interval const & rhs)
    {
        interval::merge(rhs);
        count += (dynamic_cast<ivcount const &>(rhs)).count;
    }

    bool operator==(const ivcount &other) const {
        return start == other.start && end == other.end && count == other.count;
    }

    bool operator!=(const ivcount &other) const {
        return !(*this == other);
    }

    int count;
};

std::ostream & operator<<(std::ostream & o, ivcount const & x)
{
    o << x.start << "-" << x.end << ":" << x.count;
    return o;
}

BOOST_AUTO_TEST_CASE(testIntervalList)
{
    IntervalList<ivcount> ivl;

    ivl.add(ivcount(10, 20, 1));
    ivl.add(ivcount(12, 30, 1));
    ivl.add(ivcount(32, 35, 1));
    ivl.add(ivcount(36, 37, 1));
    ivl.add(ivcount(38, 40, 1));
    ivl.add(ivcount(42, 45, 1));

    std::list<ivcount> expected = {
        ivcount(10, 30, 2),
        ivcount(32, 40, 3),
        ivcount(42, 45, 1)
    };

    listcmp(expected, ivl.getIntervals());

    BOOST_CHECK_EQUAL(ivl.query(11, 12).count, 2);
    BOOST_CHECK_EQUAL(ivl.query(31, 37).count, 3);
    BOOST_CHECK_EQUAL(ivl.query(31, 39).count, 3);
    BOOST_CHECK_EQUAL(ivl.query(42, 44).count, 1);
    BOOST_CHECK_EQUAL(ivl.query(41, 41).count, 0);
    BOOST_CHECK_EQUAL(ivl.query(45, 45).count, 1);

    ivl.keep_only(31, 44);
    BOOST_CHECK_EQUAL(ivl.query(11, 12).count, 0);
    BOOST_CHECK_EQUAL(ivl.query(31, 37).count, 3);
    BOOST_CHECK_EQUAL(ivl.query(31, 39).count, 3);
    BOOST_CHECK_EQUAL(ivl.query(42, 44).count, 1);
    BOOST_CHECK_EQUAL(ivl.query(41, 41).count, 0);
    BOOST_CHECK_EQUAL(ivl.query(45, 45).count, 0);

    expected = {
        ivcount(32, 40, 3),
        ivcount(42, 44, 1)
    };

    listcmp(expected, ivl.getIntervals());
}

BOOST_AUTO_TEST_CASE(testIntervalList2)
{
    struct ivlist : public interval
    {
        ivlist() : interval() {}
        ivlist(int64_t _start, int64_t _end) : interval(_start, _end) {
            contained_ivs.push_back(std::pair<int64_t, int64_t>(_start, _end));
        }
        virtual void merge(interval const & rhs)
        {
            contained_ivs.insert(contained_ivs.end(),
                                 (dynamic_cast<ivlist const &>(rhs)).contained_ivs.begin(),
                                 (dynamic_cast<ivlist const &>(rhs)).contained_ivs.end());
            interval::merge(rhs);
        }

        std::list< std::pair<int64_t, int64_t> > contained_ivs;
    };

    IntervalList<ivlist> ivl;

    ivl.add(ivlist(10, 20));
    ivl.add(ivlist(12, 30));
    ivl.add(ivlist(32, 35));
    ivl.add(ivlist(36, 37));
    ivl.add(ivlist(38, 40));
    ivl.add(ivlist(42, 45));

    std::list<std::pair<int64_t, int64_t> > expected, actual;

    expected = {
        std::pair<int64_t, int64_t>{10, 20},
        std::pair<int64_t, int64_t>{12, 30}
    };
    actual = ivl.query(11, 12).contained_ivs;
    listcmp(expected, actual);

    expected = {
        std::pair<int64_t, int64_t>{32, 35},
        std::pair<int64_t, int64_t>{36, 37},
        std::pair<int64_t, int64_t>{38, 40}
    };
    actual = ivl.query(31, 37).contained_ivs;
    listcmp(expected, actual);
}

BOOST_AUTO_TEST_CASE(testIntervalList3)
{
    struct ivlist : public interval
    {
        ivlist() : interval() {}
        ivlist(int64_t _start, int64_t _end) : interval(_start, _end) {
            contained_ivs.push_back(std::pair<int64_t, int64_t>(_start, _end));
        }
        virtual void merge(interval const & rhs)
        {
            contained_ivs.insert(contained_ivs.end(),
                                 (dynamic_cast<ivlist const &>(rhs)).contained_ivs.begin(),
                                 (dynamic_cast<ivlist const &>(rhs)).contained_ivs.end());
            interval::merge(rhs);
        }

        std::list< std::pair<int64_t, int64_t> > contained_ivs;
    };

    IntervalList<ivlist> ivl;

    ivl.add(ivlist(2, 4));
    ivl.add(ivlist(6, 7));

    std::list<std::pair<int64_t, int64_t> > expected, actual;

    expected = {
        std::pair<int64_t, int64_t>{2, 4},
        std::pair<int64_t, int64_t>{6, 7}
    };
    actual = ivl.query(0, 10).contained_ivs;
    listcmp(expected, actual);
}

BOOST_AUTO_TEST_CASE(testIntervalListRandom)
{
    static const int count = 2048;
    static const int icount = 20;
    static const int tcount = 300;

    for (int k = 0; k < tcount; ++k) {
        IntervalList<interval> ivl;
        bool ivs[count];
        std::fill(std::begin(ivs), std::end(ivs), false);

        for (int i = 0; i < icount; ++i) {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t) count - 1);
            std::fill(ivs + start, ivs + end + 1, true);
            ivl.add(interval(start, end));
        }

        int64_t iv_start = 0;
        int64_t iv_end = count-1;
#ifdef _DEBUG
        // remember for debugging
        IntervalList<interval> ivl_before = ivl;
#endif
        if (k % 3 == 0)
        {
            iv_end = (rand()%count) - 1;
            ivl.remove_to(iv_end);
            std::fill(ivs, ivs + iv_end + 1, false);
        }
        else if (k % 3 == 1)
        {
            iv_start = (rand()%count) - 1;
            ivl.remove_from(iv_start);
            std::fill(ivs + iv_start, std::end(ivs), false);
        }

        for (int i = 0; i < count; ++i) {
            int64_t start = rand() % count;
            int64_t end = std::min(start + rand() % 100, (int64_t) count - 1);
            bool iv_ovl = std::any_of(ivs + start, ivs + end + 1, [](bool b){return b;});
            interval iv = ivl.query(start, end);
            bool ivl_ovl = iv.start >= 0 && iv.end >= 0;
            BOOST_CHECK_EQUAL(iv_ovl, ivl_ovl);
        }
    }

}
