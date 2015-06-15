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
 * Test bulk variant preprocessing
 *
 * \file test_variantprocessing.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <sstream>

#include "Variant.hh"
#include "variant/VariantAlleleRemover.hh"
#include "variant/VariantAlleleSplitter.hh"
#include "variant/VariantLocationAggregator.hh"
#include "variant/VariantAlleleUniq.hh"

using namespace variant;

BOOST_AUTO_TEST_CASE(testProcessingAlleleRemover)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("test.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    const char * expected[] = {
        "chr1:16208-16209 0/0 PASS",
        "chr1:16210-16210 16210-16210:T 1|1 PASS",
        "chr1:16211-16211 16211-16211:T 1|1 PASS",
        "chr1:16213-16213 16213-16213:T 1|0 PASS",
        "chr1:16216-16216 16216-16216:T 1|1 PASS",
        "chr1:16217-16217 16217-16217:T 1|1 PASS",
        "chr1:16218-16218 16218-16218:A 1|1 TEST",
        "chr1:16228-16228 16228-16228:T 1|1 TEST",
        "chr1:16229-16229 16229-16229:T 1|1 PASS",
        "chr1:16255-16255 16255-16255:T 1|1 PASS",
        "chr1:16256-16256 16256-16256:T 1|1 PASS",
        "chr1:16257-16257 . PASS",  // 16257-16257:T gets removed!
        "chr1:16258-16259 0/0 PASS",
        "chr1:19994-19995 19994-19995:C 1/0 PASS",
        "chr1:19996-19997 19996-19997:C 1/1 PASS",
        "chr1:19998-19999 19998-19999:G 1/1 PASS",
        "chr1:20000-20002 20000-20002:C 1/0 PASS",
        "chr1:20004-20005 20004-20005:T 1/0 PASS",
        "chr1:826157-826160 0/0 PASS",
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS",
        "chr1:826161-826167 0/0 PASS",
        "chr1:826168-826180 0/0 PASS",
        "chr1:826180-826182 826180-826182:C 0/1 191 PASS"
    };

    VariantBufferMode bms[] = {
        VariantBufferMode::buffer_all,
        VariantBufferMode::buffer_endpos,
        VariantBufferMode::buffer_count,
        VariantBufferMode::buffer_block,
    };
    int64_t params[] = {
        0,
        16259,
        3,
        2,
    };
    size_t counts[] = {
        23,
        13,
        23,
        23,
    };

    for (int i = 0; i < 4; ++i)
    {
        VariantReader r;
        r.addSample(tp.string().c_str(), "NA12877");
        VariantProcessor pr;
        VariantAlleleRemover ar;
        pr.setReader(r, bms[i], params[i]);
        pr.addStep(ar);

        size_t count = 0;

        while(pr.advance())
        {
            Variants & v = pr.current();
            std::ostringstream ss;
            // std::cerr << v << "\n";
            ss << v;
            BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
            BOOST_CHECK_EQUAL(ss.str(), expected[count]);
            ++count;
        }
        BOOST_CHECK_EQUAL(count, counts[i]);
    }
}

BOOST_AUTO_TEST_CASE(testPreprocessingSplit)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                   .parent_path()   // src
                                    / boost::filesystem::path("example")
                                    / boost::filesystem::path("multimerge")
                                    / boost::filesystem::path("hap_alleles_");

    VariantReader r;
    r.addSample((tp.string() + "1.vcf.gz").c_str(), "SIMPLE");
    r.addSample((tp.string() + "2.vcf.gz").c_str(), "SIMPLE");

    VariantProcessor pr;
    pr.setReader(r, VariantBufferMode::buffer_block, 100);
    VariantAlleleSplitter a;
    pr.addStep(a);

    const char * expected[] = {
        "chr1:19472001-19472021 19472001-19472021:G 0/1 0/1",
        "chr1:19472001-19472021 19472001-19472021:GACACACACACACACACACACAC . 0/1",
        "chr1:19472021-19472021 19472021-19472021:CAC 0/1 .",
        "chr1:31755329-31755337 31755329-31755337:G 0/1 0/1",
        "chr1:31755329-31755337 31755329-31755337:GTATC . 0/1",
        "chr1:31755333-31755337 31755333-31755337:C 0/1 .",
    };

    int count = 0;
    while(pr.advance())
    {
        Variants & v = pr.current();
        std::ostringstream ss;
        // std::cerr << v << "\n";
        ss << v;
        BOOST_CHECK(((size_t)count) < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        ++count;
    }
    BOOST_CHECK_EQUAL(count, 6);
}

BOOST_AUTO_TEST_CASE(testPreprocessingMerge)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                   .parent_path()   // src
                                    / boost::filesystem::path("example")
                                    / boost::filesystem::path("multimerge")
                                    / boost::filesystem::path("hap_alleles_leftshifted.vcf.gz");

    VariantReader r;
    r.addSample(tp.string().c_str(), "hap_alleles_1.vcf");
    r.addSample(tp.string().c_str(), "hap_alleles_2.vcf");

    VariantProcessor pr;
    pr.setReader(r, VariantBufferMode::buffer_block, 100);
    VariantLocationAggregator a;
    a.setAggregationType(VariantLocationAggregator::aggregate_hetalt);
    pr.addStep(a);
    VariantAlleleUniq b;
    pr.addStep(b);

    const char * expected[] = {
        "chr1:19472001-19472021 19472001-19472001:GAC 19472001-19472021:G 1/2 1/2",
        "chr1:31755329-31755337 31755329-31755333:G 31755329-31755337:G 1/2 1/2",
    };

    int count = 0;
    while(pr.advance())
    {
        Variants & v = pr.current();
        std::ostringstream ss;
        // std::cerr << v << "\n";
        ss << v;
        BOOST_CHECK(((size_t)count) < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        ++count;
    }
    BOOST_CHECK_EQUAL(count, 2);
}
