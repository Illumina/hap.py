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
 * \file test_variants_cache.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <sstream>
#include <fstream>

#include "Variant.hh"
#include "helpers/StringUtil.hh"

using namespace variant;

BOOST_AUTO_TEST_CASE(testVariantStatisticsCountRefvar)
{
    VariantStatistics vs;

    RefVar rv;

    // snp
    rv.start = 10;
    rv.end = 10;
    rv.alt = "T";
    vs.add(rv);

    // del
    rv.start = 20;
    rv.end = 30;
    rv.alt = "A";
    vs.add(rv);

    // del2
    rv.start = 20;
    rv.end = 21;
    rv.alt = "A";
    vs.add(rv);

    // ins
    rv.start = 30;
    rv.end = 30;
    rv.alt = "AA";
    vs.add(rv);

    // ins
    rv.start = 40;
    rv.end = 39;
    rv.alt = "AA";
    vs.add(rv);

    // ins
    rv.start = 50;
    rv.end = 49;
    rv.alt = "A";
    vs.add(rv);

    for (int i = 0; i < 5; ++i)
    {
        // blocksubst
        rv.start = 50*i;
        rv.end = 50*i + 2;
        rv.alt = "AAA";
        vs.add(rv);
    }

    for (int i = 0; i < 10; ++i)
    {
        // complexdel
        rv.start = 50*i;
        rv.end = 50*i + 2;
        rv.alt = "AA";
        vs.add(rv);
    }

    for (int i = 0; i < 20; ++i)
    {
        // complexins
        rv.start = 50*i;
        rv.end = 50*i + 3;
        rv.alt = "AAAAAAAAAA";
        vs.add(rv);
    }

    Json::Value result(vs.write());
    BOOST_CHECK_EQUAL(result["alleles"].asUInt64(), (uint64_t)41);
    BOOST_CHECK_EQUAL(result["snp"].asUInt64(), (uint64_t)1);
    BOOST_CHECK_EQUAL(result["del"].asUInt64(), (uint64_t)2);
    BOOST_CHECK_EQUAL(result["ins"].asUInt64(), (uint64_t)3);
    BOOST_CHECK_EQUAL(result["blocksubst"].asUInt64(), (uint64_t)5);
    BOOST_CHECK_EQUAL(result["blocksubst_del"].asUInt64(), (uint64_t)10);
    BOOST_CHECK_EQUAL(result["blocksubst_ins"].asUInt64(), (uint64_t)20);
}

BOOST_AUTO_TEST_CASE(testVariantStatisticsCountVariants)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                   .parent_path()   // src
                                    / boost::filesystem::path("example") 
                                    / boost::filesystem::path("hc.vcf.gz");

    VariantReader r;
    VariantStatistics vs(true);
    int s = r.addSample(tp.c_str(), "NA12878");

    size_t count = 0;
    r.rewind();
    while(r.advance())
    {
        Variants & v = r.current();
        vs.add(v, s);
        ++count;
    }

    Json::Value result = vs.write();

    BOOST_CHECK_EQUAL(result["het"].asUInt64(), (uint64_t)2518);
    BOOST_CHECK_EQUAL(result["homalt"].asUInt64(), (uint64_t)2194);    
    BOOST_CHECK_EQUAL(result["homref"].asUInt64(), (uint64_t)79057);
    BOOST_CHECK_EQUAL(result["haploid"].asUInt64(), (uint64_t)4);
    BOOST_CHECK_EQUAL(result["hetalt"].asUInt64(), (uint64_t)66);
    BOOST_CHECK_EQUAL(result["unknown_gt"].asUInt64(), (uint64_t)246);
    BOOST_CHECK_EQUAL(result["locations"].asUInt64(), (uint64_t)84085);

    const char * vts [] = { "het", "homalt", "haploid", "hetalt" };

    for(const char * vt : vts)
    {
        uint64_t count = 0;
        for(std::string const & s : result.getMemberNames())
        {
            size_t p = s.find("__");
            if (p == std::string::npos)
            {
                continue;
            }
            if (stringutil::endsWith(s, std::string("__") + vt))
            {
                count += result[s].asUInt64();
            }
        }
        BOOST_CHECK_EQUAL(result[vt].asUInt64(), count);
    }

    BOOST_CHECK_EQUAL(result["del"].asUInt64(), (uint64_t)548);
    BOOST_CHECK_EQUAL(result["ins"].asUInt64(), (uint64_t)569);
    BOOST_CHECK_EQUAL(result["snp"].asUInt64(), (uint64_t)5888);
    BOOST_CHECK_EQUAL(result["het"].asUInt64(), (uint64_t)2518);

    BOOST_CHECK_EQUAL(result["locations"].asUInt64(), 
                      result["het"].asUInt64() + 
                      result["homalt"].asUInt64() + 
                      result["hetalt"].asUInt64() + 
                      result["homref"].asUInt64() + 
                      result["haploid"].asUInt64() + 
                      result["unknown_gt"].asUInt64());

    BOOST_CHECK(count > 0);
}
