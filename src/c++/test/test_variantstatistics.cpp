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
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path ptp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                   .parent_path()   // src
                                    / boost::filesystem::path("example")
                                    / boost::filesystem::path("chr21.fa");

    VariantStatistics vs(ptp.c_str());

    RefVar rv;

    // 1x snp
    rv.start = 10;
    rv.end = 10;
    rv.alt = "T";
    vs.add("chr21", rv);

    // 10x del
    rv.start = 20;
    rv.end = 30;
    rv.alt = "N";
    vs.add("chr21", rv);

    // 1x SNP - 1x del (complex_sd)
    rv.start = 20;
    rv.end = 21;
    rv.alt = "A";
    vs.add("chr21", rv);

    // 1x SNP - 1x ins (complex_si)
    rv.start = 30;
    rv.end = 30;
    rv.alt = "AA";
    vs.add("chr21", rv);

    // 2x ins
    rv.start = 40;
    rv.end = 39;
    rv.alt = "AA";
    vs.add("chr21", rv);

    // 1x ins
    rv.start = 50;
    rv.end = 49;
    rv.alt = "A";
    vs.add("chr21", rv);

    for (int i = 0; i < 5; ++i)
    {
        // blocksubst -- 3xSNP x 5
        rv.start = 50*i;
        rv.end = 50*i + 2;
        rv.alt = "AAA";
        vs.add("chr21", rv);
    }

    for (int i = 0; i < 10; ++i)
    {
        // complexdel -- 2x SNP 1x del x10
        rv.start = 50*i;
        rv.end = 50*i + 2;
        rv.alt = "AA";
        vs.add("chr21", rv);
    }

    for (int i = 0; i < 20; ++i)
    {
        // complexins -- 4x snp, 6 x ins x20
        rv.start = 50*i;
        rv.end = 50*i + 3;
        rv.alt = "AAAAAAAAAA";
        vs.add("chr21", rv);
    }

    Json::Value result(vs.write());

//    Json::FastWriter w;
//    std::cerr << w.write(result) << std::endl;

    BOOST_CHECK_EQUAL(result["alleles__snp"].asUInt64(), (uint64_t)6); // 5 + 1 pure SNP records
    BOOST_CHECK_EQUAL(result["alleles__ins"].asUInt64(), (uint64_t)2); // 2 pure insertion records
    BOOST_CHECK_EQUAL(result["alleles__del"].asUInt64(), (uint64_t)1); // 1 pure deletion
    BOOST_CHECK_EQUAL(result["alleles__complex_sd"].asUInt64(), (uint64_t)(1 + 10));
    BOOST_CHECK_EQUAL(result["alleles__complex_si"].asUInt64(), (uint64_t)(1 + 20));
    BOOST_CHECK_EQUAL(result["nucleotides__snp"].asUInt64(), (uint64_t)(1 + 1 + 1 + 5*3 + 2*10 + 4*20));
    BOOST_CHECK_EQUAL(result["nucleotides__ins"].asUInt64(), (uint64_t)(1 + 2 + 1 + 6*20));
    BOOST_CHECK_EQUAL(result["nucleotides__del"].asUInt64(), (uint64_t)(10 + 1 + 10));
}

BOOST_AUTO_TEST_CASE(testVariantStatisticsCountVariants)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path ptp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                   .parent_path()   // src
                                    / boost::filesystem::path("example");
    boost::filesystem::path tp = ptp / boost::filesystem::path("hc.vcf.gz");
    boost::filesystem::path ftp = ptp / boost::filesystem::path("chr21.fa");

    VariantReader r;
    VariantStatistics vs(ftp.c_str(), true);
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
    Json::FastWriter w;
    std::cerr << w.write(result) << std::endl;

    //{"alleles__del":426,
    // "alleles__ins":431,
    // "alleles__snp":3991,
    // "hemi__snp":4,
    // "het__nocall":2088,
    // "het__unknown":223,
    // "hetalt__complex_id":12,
    // "hetalt__del":169,
    // "hetalt__ins":181,
    // "hetalt__snp":1898,
    // "hom__ref":79057,
    // "nocall__none":246,
    // "nucleotides__del":1329,
    // "nucleotides__ins":1169,
    // "nucleotides__ref":1022,
    // "nucleotides__snp":3991}

    BOOST_CHECK_EQUAL(result["nucleotides__del"].asUInt64(), (uint64_t)1329);
    BOOST_CHECK_EQUAL(result["nucleotides__ins"].asUInt64(), (uint64_t)1169);
    BOOST_CHECK_EQUAL(result["nucleotides__snp"].asUInt64(), (uint64_t)3991);
    BOOST_CHECK_EQUAL(result["nucleotides__ref"].asUInt64(), (uint64_t)1022);

    // TODO put these back in after validating by outputting with VCF
//    BOOST_CHECK_EQUAL(result["het"].asUInt64(), (uint64_t)2518);
//    BOOST_CHECK_EQUAL(result["homalt"].asUInt64(), (uint64_t)2194);
//    BOOST_CHECK_EQUAL(result["homref"].asUInt64(), (uint64_t)79057);
//    BOOST_CHECK_EQUAL(result["haploid"].asUInt64(), (uint64_t)4);
//    BOOST_CHECK_EQUAL(result["hetalt"].asUInt64(), (uint64_t)66);
//    BOOST_CHECK_EQUAL(result["unknown_gt"].asUInt64(), (uint64_t)246);
//    BOOST_CHECK_EQUAL(result["locations"].asUInt64(), (uint64_t)84085);

    BOOST_CHECK(count > 0);
}
