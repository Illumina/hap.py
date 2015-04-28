// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Adapted from interval_tree_test.cpp from https://github.com/ekg/intervaltree
// Converted to boost test framework.
// 
// Copyright (c) 2011 Erik Garrison

// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/**
 *  \brief Interval tree test cases
 *
 *
 * \file test_haplotypes.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include <iostream>
#include <sstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <random>
#include <time.h>
#include <assert.h>
#include "IntervalTree.h"

using namespace std;

typedef Interval<bool> interval;
typedef vector<interval> intervalVector;
typedef IntervalTree<bool> intervalTree;
typedef vector<std::size_t> countsVector;

template<typename K>
K randKey(K floor, K ceiling) {
    K range = ceiling - floor;
    return floor + range * ((double) rand() / (double) (RAND_MAX + 1.0));
}

template<class T, typename K>
Interval<T,K> randomInterval(K maxStart, K maxLength, K maxStop, const T& value) {
    K start = randKey<K>(0, maxStart);
    K stop = min<K>(randKey<K>(start, start + maxLength), maxStop);
    return Interval<T,K>(start, stop, value);
}

BOOST_AUTO_TEST_CASE(testIntervalTreeSanitySimple)
{
    // a simple sanity check
    intervalVector sanityIntervals;
    sanityIntervals.push_back(interval(60, 80, true));
    sanityIntervals.push_back(interval(20, 40, true));
    intervalTree sanityTree(sanityIntervals);

    intervalVector sanityResults;
    sanityTree.findOverlapping(30, 50, sanityResults);
    BOOST_CHECK_EQUAL(sanityResults.size(), (size_t)1);

    sanityResults.clear();
    sanityTree.findContained(15, 45, sanityResults);
    BOOST_CHECK_EQUAL(sanityResults.size(), (size_t)1);
}

BOOST_AUTO_TEST_CASE(testIntervalTreeSpeed)
{
    srand((unsigned)time(NULL));

    intervalVector intervals;
    intervalVector queries;
    
    // generate a test set of target intervals
    for (int i = 0; i < 10000; ++i) 
    {
        intervals.push_back(randomInterval<bool>(100000, 1000, 100000 + 1, true));
    }
    
    // and queries
    for (int i = 0; i < 5000; ++i) 
    {
        queries.push_back(randomInterval<bool>(100000, 1000, 100000 + 1, true));
    }

    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;

    // using brute-force search
    countsVector bruteforcecounts;
    Clock::time_point t0 = Clock::now();
    for (intervalVector::iterator q = queries.begin(); q != queries.end(); ++q) 
    {
        intervalVector results;
        for (intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) 
        {
            if (i->start >= q->start && i->stop <= q->stop) 
            {
                results.push_back(*i);
            }
        }
        bruteforcecounts.push_back(results.size());
    }

    Clock::time_point t1 = Clock::now();
    milliseconds ms = chrono::duration_cast<milliseconds>(t1 - t0);
    cout << "brute force:\t" << ms.count() << "ms" << endl;

    // using the interval tree
    intervalTree tree = intervalTree(intervals);
    countsVector treecounts;
    t0 = Clock::now();
    for (intervalVector::iterator q = queries.begin(); q != queries.end(); ++q) 
    {
        intervalVector results;
        tree.findContained(q->start, q->stop, results);
        treecounts.push_back(results.size());
    }
    t1 = Clock::now();
    ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
    cout << "interval tree:\t" << ms.count() << "ms" << endl;

    // check that the same number of results are returned
    countsVector::iterator b = bruteforcecounts.begin();
    for (countsVector::iterator t = treecounts.begin(); t != treecounts.end(); ++t, ++b) 
    {
        BOOST_CHECK_EQUAL(*b, *t);
    }
}

