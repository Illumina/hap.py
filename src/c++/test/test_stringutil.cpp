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
 * \file test_stringutil.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "helpers/StringUtil.hh"

BOOST_AUTO_TEST_CASE(testStringSplit)
{
	std::vector<std::string> v;

	stringutil::split("a, b, c", v);

    BOOST_CHECK(v.size() == 3);
    BOOST_CHECK(v[0] == "a");
    BOOST_CHECK(v[1] == "b");
    BOOST_CHECK(v[2] == "c");

    v.clear();

	stringutil::split("a, b, c", v, " ,", true);

    BOOST_CHECK(v.size() == 5);
    BOOST_CHECK(v[0] == "a");
    BOOST_CHECK(v[1] == "");
    BOOST_CHECK(v[2] == "b");
    BOOST_CHECK(v[3] == "");
    BOOST_CHECK(v[4] == "c");

    v.clear();

    stringutil::split("a: b- c", v, " :-");

    BOOST_CHECK(v.size() == 3);
    BOOST_CHECK(v[0] == "a");
    BOOST_CHECK(v[1] == "b");
    BOOST_CHECK(v[2] == "c");
}

BOOST_AUTO_TEST_CASE(testStringEndsWith)
{
    BOOST_CHECK(stringutil::endsWith("abc", "bc"));
    BOOST_CHECK(!stringutil::endsWith("abc", "cb"));
}

BOOST_AUTO_TEST_CASE(testStringReplaceAll)
{
    BOOST_CHECK_EQUAL(stringutil::replaceAll("1,000,000", ",", ""), "1000000");
    BOOST_CHECK_EQUAL(stringutil::replaceAll("1,000,000", ",0", "2"), "1200200");
}

BOOST_AUTO_TEST_CASE(testStringFormatPos)
{
    BOOST_CHECK_EQUAL(stringutil::formatPos("chr1"), "chr1");
    BOOST_CHECK_EQUAL(stringutil::formatPos("chr1", 999), "chr1:1000");
    BOOST_CHECK_EQUAL(stringutil::formatPos("chr1", 999, 1999), "chr1:1000-2000");
}

BOOST_AUTO_TEST_CASE(testStringParsePos)
{
    std::string chr = "-";
    int64_t start = -1, end = -1;

    stringutil::parsePos("chr1", chr, start, end);

    BOOST_CHECK_EQUAL(chr, "chr1");
    BOOST_CHECK_EQUAL(start, -1);
    BOOST_CHECK_EQUAL(end, -1);

    chr = "-";
    start = -1;
    end = -1;

    stringutil::parsePos("chr1:1,000", chr, start, end);

    BOOST_CHECK_EQUAL(chr, "chr1");
    BOOST_CHECK_EQUAL(start, 999);
    BOOST_CHECK_EQUAL(end, -1);

    chr = "-";
    start = -1;
    end = -1;

    stringutil::parsePos("chr1:1,000-2000", chr, start, end);

    BOOST_CHECK_EQUAL(chr, "chr1");
    BOOST_CHECK_EQUAL(start, 999);
    BOOST_CHECK_EQUAL(end, 1999);
}
