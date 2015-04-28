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
 * Test cases for diploid reference class
 *
 * \file test_diploidreference.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/filesystem/path.hpp>

#include "DiploidReference.hh"

#include <iostream>
#include <sstream>
#include <vector>
#include <set>

using namespace variant;
using namespace haplotypes;

BOOST_AUTO_TEST_CASE(diploidReferenceBasic)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data");

    std::string datapath = tp.string();

    DiploidReference dr((datapath + "/refgraph1.vcf.gz").c_str(), 
                        "NA12877", (datapath + "/chrQ.fa").c_str());

    dr.setRegion("chrQ", 0, 25);

    std::set<std::string> expected;

    // A[AAC/]CC[GGG]AAA[C/G] [C|T][C|T][T|G]AACCCGG[C|T]TTTG
    // -> haploid 
    // ACCGGGAAACCCTAACCCGGCTTTGG
    // ACCGGGAAACTTGAACCCGGTTTTGG
    // ACCGGGAAAGCCTAACCCGGCTTTGG
    // ACCGGGAAAGTTGAACCCGGTTTTGG
    // AAACCCGGGAAACCCTAACCCGGCTTTGG
    // AAACCCGGGAAACTTGAACCCGGTTTTGG
    // AAACCCGGGAAAGCCTAACCCGGCTTTGG
    // AAACCCGGGAAAGTTGAACCCGGTTTTGG

    expected.insert("hetalt(ACCGGGAAACCCTAACCCGGCTTTGG|AAACCCGGGAAAGTTGAACCCGGTTTTGG)");
    expected.insert("hetalt(ACCGGGAAACTTGAACCCGGTTTTGG|AAACCCGGGAAAGCCTAACCCGGCTTTGG)");
    expected.insert("hetalt(ACCGGGAAAGCCTAACCCGGCTTTGG|AAACCCGGGAAACTTGAACCCGGTTTTGG)");
    expected.insert("hetalt(ACCGGGAAAGTTGAACCCGGTTTTGG|AAACCCGGGAAACCCTAACCCGGCTTTGG)");

    size_t count = 0;
    while(dr.hasNext())
    {
        std::ostringstream ss;
        ss << dr.next();
        // std::cerr << ss.str() << "\n";
        BOOST_CHECK_EQUAL(expected.count(ss.str()), (size_t)1);
        dr.advance();
        ++count;
    }
    BOOST_CHECK_EQUAL(count, expected.size());
}

BOOST_AUTO_TEST_CASE(diploidReferenceBasicHet)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data");

    std::string datapath = tp.string();

    DiploidReference dr((datapath + "/refgraph1.vcf.gz").c_str(), 
                        "NA12877", (datapath + "/chrQ.fa").c_str());
        

    dr.setRegion("chrQ", 0, 4);

    std::set<std::string> expected;

    expected.insert("het(AC|AAACC)");

    size_t count = 0;
    while(dr.hasNext())
    {
        std::ostringstream ss;
        ss << dr.next();
        BOOST_CHECK_EQUAL(expected.count(ss.str()), (size_t)1);
        dr.advance();
        ++count;
    }
    BOOST_CHECK_EQUAL(count, expected.size());
}

BOOST_AUTO_TEST_CASE(diploidReferencePhased)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data");

    std::string datapath = tp.string();

    DiploidReference dr((datapath + "/refgraph1.vcf.gz").c_str(), 
                        "NA12877", (datapath + "/chrQ.fa").c_str());
        

    dr.setRegion("chrQ", 11, 25);

    size_t count = 0;
    while(dr.hasNext())
    {
        std::ostringstream ss;
        ss << dr.next();
        // accgggaaagtTGAACCCGGTTTTGG|aaacccgggaaaccCTAACCCGGCTTTGG
        BOOST_CHECK_EQUAL(ss.str(), "hetalt(TGAACCCGGTTTTGG|CTAACCCGGCTTTGG)");
        dr.advance();
        ++count;
    }
    BOOST_CHECK_EQUAL(count, (size_t)1);
}

BOOST_AUTO_TEST_CASE(diploidReferenceHomref)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data");

    std::string datapath = tp.string();

    DiploidReference dr((datapath + "/refgraph1.vcf.gz").c_str(), 
                        "NA12877", (datapath + "/chrQ.fa").c_str());

    // region without variation should return homref
    dr.setRegion("chrQ", 27, 36);
    size_t count = 0;
    while(dr.hasNext())
    {
        std::ostringstream ss;
        ss << dr.next();
        BOOST_CHECK_EQUAL(ss.str(), "homref(TTTGGGTTT)");
        dr.advance();
        ++count;
    }
    BOOST_CHECK_EQUAL(count, (size_t)1);
}