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

    FastaFile f(ptp.c_str());
    VariantStatistics vs(f);

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

    BOOST_CHECK_EQUAL(result["al__s"].asUInt64(), (uint64_t)6); // 5 + 1 pure SNP records
    BOOST_CHECK_EQUAL(result["al__i"].asUInt64(), (uint64_t)2); // 2 pure insertion records
    BOOST_CHECK_EQUAL(result["al__d"].asUInt64(), (uint64_t)1); // 1 pure deletion
    BOOST_CHECK_EQUAL(result["al__sd"].asUInt64(), (uint64_t)(1 + 10));
    BOOST_CHECK_EQUAL(result["al__si"].asUInt64(), (uint64_t)(1 + 20));
    BOOST_CHECK_EQUAL(result["nuc__s"].asUInt64(), (uint64_t)(1 + 1 + 1 + 5*3 + 2*10 + 4*20));
    BOOST_CHECK_EQUAL(result["nuc__i"].asUInt64(), (uint64_t)(1 + 2 + 1 + 6*20));
    BOOST_CHECK_EQUAL(result["nuc__d"].asUInt64(), (uint64_t)(10 + 1 + 10));
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
    FastaFile f(ftp.c_str());
    VariantStatistics vs(f, true);
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
    /* bcftools stats ../hap.py/example/hc.vcf.gz | grep ^SN */
    /* SN      0       number of samples:      1 */
    /* SN      0       number of records:      84085 */
    /* SN      0       number of SNPs: 3990 */
    /* SN      0       number of MNPs: 0 */
    /* SN      0       number of indels:       792 */
    /* SN      0       number of others:       0 */
    /* SN      0       number of multiallelic sites:   66 */
    /* SN      0       number of multiallelic SNP sites:       1 */
    /* Json::FastWriter w; */
    /* std::cerr << w.write(result) << std::endl; */
    //{"al__d":426,  426 + 431 - 65 het-alt locations = 792 indels
    // "al__i":431,
    // "al__s":3991, 3991 - 1 hetalt site = 3990
    // "hemi__s":4,
    // "het__rd":223,
    // "het__ri":207,
    // "het__rs":2088,
    // "hetalt__d":22,  # of multiallelic sites: 22+31+12+1 = 66
    // "hetalt__i":31,
    // "hetalt__id":12,
    // "hetalt__s":1,
    // "homalt__d":147,
    // "homalt__i":150,
    // "homalt__s":1897,
    // "homref__r":79057,
    // "nocall__nc":246,
    // "nuc__d":1329,
    // "nuc__i":1169,
    // "nuc__r":1022,
    // "nuc__s":3991}
    BOOST_CHECK_EQUAL(result["nuc__d"].asUInt64(), (uint64_t)1329);
    BOOST_CHECK_EQUAL(result["nuc__i"].asUInt64(), (uint64_t)1169);
    BOOST_CHECK_EQUAL(result["nuc__s"].asUInt64(), (uint64_t)3991);
    BOOST_CHECK_EQUAL(result["nuc__r"].asUInt64(), (uint64_t)1022);

    BOOST_CHECK_EQUAL(result["al__d"].asUInt64(), (uint64_t)426);
    BOOST_CHECK_EQUAL(result["al__i"].asUInt64(), (uint64_t)431);
    BOOST_CHECK_EQUAL(result["al__s"].asUInt64(), (uint64_t)3991);
    BOOST_CHECK_EQUAL(result["hemi__s"].asUInt64(), (uint64_t)4);
    BOOST_CHECK_EQUAL(result["het__rd"].asUInt64(), (uint64_t)223);
    BOOST_CHECK_EQUAL(result["het__ri"].asUInt64(), (uint64_t)207);
    BOOST_CHECK_EQUAL(result["het__rs"].asUInt64(), (uint64_t)2088);
    BOOST_CHECK_EQUAL(result["hetalt__d"].asUInt64(), (uint64_t)22);
    BOOST_CHECK_EQUAL(result["hetalt__i"].asUInt64(), (uint64_t)31);
    BOOST_CHECK_EQUAL(result["hetalt__id"].asUInt64(), (uint64_t)12);
    BOOST_CHECK_EQUAL(result["hetalt__s"].asUInt64(), (uint64_t)1);
    BOOST_CHECK_EQUAL(result["homalt__d"].asUInt64(), (uint64_t)147);
    BOOST_CHECK_EQUAL(result["homalt__i"].asUInt64(), (uint64_t)150);
    BOOST_CHECK_EQUAL(result["homalt__s"].asUInt64(), (uint64_t)1897);
    BOOST_CHECK_EQUAL(result["homref__r"].asUInt64(), (uint64_t)79057);
    BOOST_CHECK_EQUAL(result["nocall__nc"].asUInt64(), (uint64_t)246);

    BOOST_CHECK(count > 0);
}
