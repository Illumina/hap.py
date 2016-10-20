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
 * \brief
 *
 * \file test_format_extraction.cpp
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
#include "helpers/BCFHelpers.hh"

using namespace variant;


BOOST_AUTO_TEST_CASE(variantReadingFormats)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("test_formats.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");

    size_t count = 0;

    const char * expected[] = {
        "chr1:826157-826160 826157-826160:C 0/0 PASS",
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS",
        "chr1:826161-826167 826161-826167:T 0/0 PASS",
        "chr1:826168-826180 826168-826180:T 0/0 PASS",
        "chr1:826180-826182 826180-826182:C 826180-826182:T 0/1 191 PASS"
    };

    float expected_qual[] = {
    	-1,
    	153.0,
    	-1,
    	-1,
		191,
    };

    float expected_gq[] = {
    	138.0,
		193.0,
		102.0,
		96.0,
		231.0,
    };

    float expected_ad[] = {
    	0, 0,
		34, 5,
    	0, 0,
    	0, 0,
		38, 7,
    };

    float expected_ad_r[] = {
    	-1,
	 	34,
    	-1,
    	-1,
		 38,
    };

    float expected_ad_o[] = {
    	0,
	 	0,
    	0,
    	0,
		0,
    };


    int expected_dp[] = {
    	47,
		55,
		35,
		33,
		48,
    };

    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream ss;
        ss << v;
        // std::cerr << ss.str() << "\n";
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);

        if(expected_qual[count] >= 0)
        {
        	BOOST_CHECK_CLOSE(v.calls[0].qual, expected_qual[count], 0.01);
        }
        else
        {
        	BOOST_CHECK(std::isnan(v.calls[0].qual));
        }

        float gq = v.calls[0].formats["GQX"].asFloat();
		if(!gq) gq = v.calls[0].formats["GQ"].asFloat();
        BOOST_CHECK_CLOSE(gq, expected_gq[count], 0.01);
        BOOST_CHECK_EQUAL(v.calls[0].ad[0], expected_ad[2*count]);
        BOOST_CHECK_EQUAL(v.calls[0].ad[1], expected_ad[2*count+1]);
        BOOST_CHECK_EQUAL(v.calls[0].ad_ref, expected_ad_r[count]);
        BOOST_CHECK_EQUAL(v.calls[0].ad_other, expected_ad_o[count]);
        BOOST_CHECK_EQUAL(v.calls[0].dp, expected_dp[count]);

        ++count;
    }
    BOOST_CHECK(count > 0);
}

BOOST_AUTO_TEST_CASE(variantReadingFormats2)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("test_formats.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");
    r.addSample(tp.string().c_str(), "NA12878");

    size_t count = 0;

    const char * expected[] = {
        "chr1:826157-826160 826157-826160:C 0/0 PASS 1/0 PASS",
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS 1/0 153 PASS",
        "chr1:826161-826167 826161-826167:T 0/0 PASS 1/0 PASS",
        "chr1:826168-826180 826168-826180:T 0/0 PASS 1/0 PASS",
        "chr1:826180-826182 826180-826182:C 826180-826182:T 0/1 191 PASS 1/0 191 PASS"
    };

    float expected_qual[] = {
    	-1,
    	153.0,
    	-1,
    	-1,
		191,
    	-1,
    	153.0,
    	-1,
    	-1,
		191,
    };

    float expected_gq[] = {
    	138.0,
		193.0,
		102.0,
		96.0,
		231.0,
    	139.0,
		194.0,
		103.0,
		97.0,
		232.0,
    };

    float expected_ad[] = {
    	0, 0,
		34, 5,
    	0, 0,
    	0, 0,
		38, 7,
    	0, 0,
		5, 34,		// GT is the other way around in this one
    	0, 0,
    	0, 0,
		8, 39,
    };

    float expected_ad_r[] = {
    	-1,
	 	34,
    	-1,
    	-1,
		 38,
    	-1,
	 	34,
    	-1,
    	-1,
		 39,
    };

    float expected_ad_o[] = {
    	0,
	 	0,
    	0,
    	0,
		0,
    	0,
	 	0,
    	0,
    	0,
		1,
    };

    int expected_dp[] = {
    	47,
		55,
		35,
		33,
		48,
    	47,
		55,
		36,
		37,
		49,
    };

    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream ss;
        ss << v;
        // std::cerr << ss.str() << "\n";
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);

        for (int i = 0; i < 2; ++i)
        {
        	// std::cerr << "\t" << i << "\n";
	        if(expected_qual[count] >= 0)
	        {
	        	BOOST_CHECK_CLOSE(v.calls[i].qual, expected_qual[count + 5*i], 0.01);
	        }
	        else
	        {
	        	BOOST_CHECK(std::isnan(v.calls[i].qual));
	        }
            float gq = v.calls[i].formats["GQX"].asFloat();
            if(gq <= 0) gq = v.calls[i].formats["GQ"].asFloat();
	        BOOST_CHECK_CLOSE(gq, expected_gq[count + 5*i], 0.01);
	        BOOST_CHECK_EQUAL(v.calls[i].ad[0], expected_ad[2*count + 10*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].ad[1], expected_ad[2*count+1 + 10*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].ad_ref, expected_ad_r[count + 5*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].ad_other, expected_ad_o[count + 5*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].dp, expected_dp[count + 5*i]);
        }
        ++count;
    }
    BOOST_CHECK(count > 0);
}
