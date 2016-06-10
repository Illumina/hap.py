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
 * \file test_variants.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <sstream>
#include <cmath>

#include "Variant.hh"

using namespace variant;

BOOST_AUTO_TEST_CASE(variantGTT)
{
    Call var;

    var.ngt = 1;
    var.gt[0] = -1;
    BOOST_CHECK(getGTType(var) == gt_unknown);
    var.gt[0] = 0;
    BOOST_CHECK(getGTType(var) == gt_homref);
    var.gt[0] = 1;
    BOOST_CHECK(getGTType(var) == gt_haploid);

    var.ngt = 2;
    var.gt[0] = 0;
    var.gt[1] = -1;
    BOOST_CHECK(getGTType(var) == gt_unknown);
    var.gt[0] = 0;
    var.gt[1] = 0;
    BOOST_CHECK(getGTType(var) == gt_homref);
    var.gt[0] = 1;
    var.gt[1] = 0;
    BOOST_CHECK(getGTType(var) == gt_het);
    var.gt[0] = 0;
    var.gt[1] = 1;
    BOOST_CHECK(getGTType(var) == gt_het);
    var.gt[0] = 1;
    var.gt[1] = 1;
    BOOST_CHECK(getGTType(var) == gt_homalt);
    var.gt[0] = 2;
    var.gt[1] = 2;
    BOOST_CHECK(getGTType(var) == gt_homalt);
    var.gt[0] = 1;
    var.gt[1] = 2;
    BOOST_CHECK(getGTType(var) == gt_hetalt);

    var.ngt = 3;
    BOOST_CHECK(getGTType(var) == gt_unknown);
}

BOOST_AUTO_TEST_CASE(variantInfo)
{
    // TODO re-add this test-case
}
 
BOOST_AUTO_TEST_CASE(variantReading)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("test.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");

    size_t count = 0;

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
        "chr1:16257-16257 16257-16257:T . PASS",
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

    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream ss;
        ss << v;
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        ++count;
    }

}


BOOST_AUTO_TEST_CASE(variantReadingRegion)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("test.vcf.gz");

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");

    r.rewind("chr1", 16257);

    size_t count = 0;

    const char * expected[] = {
        "chr1:16257-16257 16257-16257:T . PASS",
        "chr1:16258-16259 0/0 PASS",
        "chr1:19994-19995 19994-19995:C 1/0 PASS",
        "chr1:19996-19997 19996-19997:C 1/1 PASS",
        "chr1:19998-19999 19998-19999:G 1/1 PASS",
        "chr1:20000-20002 20000-20002:C 1/0 PASS",
    };

    while(r.advance())
    {
        Variants & v = r.current();

        BOOST_CHECK_EQUAL(v.chr, "chr1");
        BOOST_CHECK(v.pos >= 16257);
        if(v.pos > 20000)
        {
            break;
        }

        std::ostringstream ss;
        ss << v;
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        ++count;
    }

}

BOOST_AUTO_TEST_CASE(variantReading2s)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("hapblocks.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");
    r.addSample(tp.string().c_str(), "NA12878");

    size_t count = 0;

    const char * expected[] = {
        "chr1:8999-8999 8999-8999:A 1/1 PASS 0/1 PASS",
        "chr1:9013-9013 9013-9013:A 1/1 PASS 0/1 PASS",
        "chr1:9499-9499 9499-9499:A 1/1 PASS 0/1 PASS",
        "chr1:9999-10002 9999-10002:A 9999-10002:AAA 0/1 PASS 0/2 PASS",
        "chr1:10004-10005 10004-10005:AGGGC 10004-10005:A 0/1 PASS 0/2 PASS",
        "chr1:10009-10009 10009-10009:G 10009-10009:T 0/1 PASS 0/2 PASS",
    };

    while(r.advance())
    {
        Variants & v = r.current();

        std::ostringstream ss;
        ss << v;
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        ++count;
    }

}

BOOST_AUTO_TEST_CASE(variantReadingRegions)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("test.vcf.gz");
    boost::filesystem::path tp2 = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("testr.bed");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.setRegions(tp2.c_str(), true);
    r.addSample(tp.c_str(), "NA12877");

    size_t count = 0;

    const char * expected[] = {
        "chr1:16211-16211 16211-16211:T 1|1 PASS",
        "chr1:16213-16213 16213-16213:T 1|0 PASS",
        "chr1:16216-16216 16216-16216:T 1|1 PASS",
        "chr1:16217-16217 16217-16217:T 1|1 PASS",
        "chr1:826157-826160 0/0 PASS",
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS",
        "chr1:826161-826167 0/0 PASS",
    };

    const std::function<void(Variants &)> expected_info[] = {
        [](Variants & ) {},
        [](Variants & ) {},
        [](Variants & ) {},
        [](Variants & ) {},
        [](Variants & v) { BOOST_CHECK_EQUAL(v.getInfoInt("END"), 826161); },
        [](Variants & v) {
            BOOST_CHECK_EQUAL(v.getInfoInt("REFREP"), 0);
            BOOST_CHECK_EQUAL(v.getInfoInt("IDREP"), 1);
            BOOST_CHECK_EQUAL(v.getInfoString("CIGAR"), "1M6I");
            BOOST_CHECK_EQUAL(v.getInfoString("RU"), "TTTGAT");
        },
        [](Variants & v) { BOOST_CHECK_EQUAL(v.getInfoInt("END"), 826168); },
    };

    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream ss;
        ss << v;
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        expected_info[count](v);
        ++count;
    }
}

BOOST_AUTO_TEST_CASE(variantReadingTargets)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("test.vcf.gz");
    boost::filesystem::path tp2 = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("testr.bed");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.setTargets(tp2.c_str(), true);
    r.addSample(tp.c_str(), "NA12877");

    size_t count = 0;

    const char * expected[] = {
        "chr1:16211-16211 16211-16211:T 1|1 PASS",
        "chr1:16213-16213 16213-16213:T 1|0 PASS",
        "chr1:16216-16216 16216-16216:T 1|1 PASS",
        "chr1:16217-16217 16217-16217:T 1|1 PASS",
        // "chr1:826157-826160 0/0 PASS", not found by targets
        // 
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS",
        "chr1:826161-826167 0/0 PASS",
    };

    while(r.advance())
    {
        Variants & v = r.current();
        // std::cerr << v << "\n";
        std::ostringstream ss;
        ss << v;
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);
        ++count;
    }

}


BOOST_AUTO_TEST_CASE(variantReadingMultisample)
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
    std::cerr << "Reading " << (tp.string() + "1.vcf.gz").c_str() << "\n";
    std::cerr << "Reading " << (tp.string() + "2.vcf.gz").c_str() << "\n";
    r.addSample((tp.string() + "1.vcf.gz").c_str(), "SIMPLE");
    r.addSample((tp.string() + "2.vcf.gz").c_str(), "SIMPLE");

    const char * expected[] = {
        "chr1:19472001-19472021 19472001-19472021:G 0/1 .",
        "chr1:19472001-19472021 19472001-19472021:GACACACACACACACACACACAC 19472001-19472021:G . 1/2",
        "chr1:19472021-19472021 19472021-19472021:CAC 0/1 .",
        "chr1:31755329-31755337 31755329-31755337:G 0/1 .",
        "chr1:31755329-31755337 31755329-31755337:GTATC 31755329-31755337:G . 1/2",
        "chr1:31755333-31755337 31755333-31755337:C 0/1 .",
    };

    int count = 0;
    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream oss;
        oss << v;
        // std::cerr << v << "\n";
        BOOST_CHECK(count < 6);
        BOOST_CHECK_EQUAL(expected[count], oss.str());
        ++count;
    }
    BOOST_CHECK(count == 6);
}

