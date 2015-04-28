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
 *  \brief Test cases for haplotype creation
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

#include <boost/filesystem/path.hpp>

#include <iostream>
#include <sstream>

#include "Haplotype.hh"

using namespace variant;
using namespace haplotypes;

BOOST_AUTO_TEST_CASE(haplotypeBasic)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrQ", tp.string().c_str());

    h.addVar(5, 9, "AC");
    BOOST_CHECK_EQUAL(h.seq(), "AC");
    int64_t min=-1, max=-1, sum=-1;
    h.lengthDiffs(min, max, sum);
    BOOST_CHECK_EQUAL(min, -3);
    BOOST_CHECK_EQUAL(max, -3);
    BOOST_CHECK_EQUAL(sum, -3);
}

BOOST_AUTO_TEST_CASE(haplotypeDel)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrQ", tp.string().c_str());

    h.addVar(5, 9, "AC");
    h.addVar(12, 13, "T");

    // >chrQ
    // AAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTT
    //      **     *
    // AAACCAC---CCT-ACCCGGGTTTGGGTTTGGGTTT

    BOOST_CHECK_EQUAL(h.seq(), "ACCCT");
    int64_t min=-1, max=0, sum=-1;
    h.lengthDiffs(min, max, sum);
    BOOST_CHECK_EQUAL(min, -4);
    BOOST_CHECK_EQUAL(max, -3);
    BOOST_CHECK_EQUAL(sum, -4);
}

BOOST_AUTO_TEST_CASE(haplotypeIns)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrQ", tp.string().c_str());

    h.addVar(5, 5, "CGCT");
    h.addVar(10, 10, "CA");

    // >chrQ
    // AAACCC---AAACC-CAAACCCGGGTTTGGGTTTGGGTTT
    //      |     
    // AAACCCGCTAAACCAC...

    BOOST_CHECK_EQUAL(h.seq(), "CGCTAAACCA");
    int64_t min=-1, max=-1, sum=-1;
    h.lengthDiffs(min, max, sum);
    BOOST_CHECK_EQUAL(min, 3);
    BOOST_CHECK_EQUAL(max, 4);
    BOOST_CHECK_EQUAL(sum, 4);
}

BOOST_AUTO_TEST_CASE(haplotypeIns2)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrQ", tp.string().c_str());

    // test inserting with zero ref length
    h.addVar(6, 5, "GCT");
    h.addVar(11, 10, "A");

    // >chrQ
    // AAACCC---AAACC-CAAACCCGGGTTTGGGTTTGGGTTT
    //      |     
    // AAACCCGCTAAACCAC...

    BOOST_CHECK_EQUAL(h.seq(5, 15), "CGCTAAACCACAAAC");

    int64_t min=-1, max=-1, sum=-1;
    h.lengthDiffs(min, max, sum);
    BOOST_CHECK_EQUAL(min, 3);
    BOOST_CHECK_EQUAL(max, 4);
    BOOST_CHECK_EQUAL(sum, 4);
}


BOOST_AUTO_TEST_CASE(haplotypeIns3)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrS", tp.string().c_str());

    // DEBUG REF: TAATGACAGCGACTTGAGACATACACCTGCGCCTAATCACTTCAGAGGTGTTGTAGACGGGAAATAAGGATCTACCCTTACTTTGTGTCACTACTAGTGAAGCTGGCTCATGTGGGGACTAGGAAGAGCTGCCTCTCGTGGGAGCGATGC
    // DEBUG ALT: TAATGACAGCGACTTGAGACATACACCTGCGCCTAATCACTTCAGAGGTGTATTGTCAGTCTGGTTAGAGCCACATGGCCTACCAGCTGACCAGTCCTTGCTTGTAGACGGGAAATAAGGATCTACCCTTACTTTGTGTCACTACTAGTGAAGCTGGCTCATGTGGGGACTAGGAAGAGCTGCCTCTCGTGGGAGCGATGC
    // DEBUG ALN: 0 149 0 200 : 51M51D99M
    // DEBUG RV: 51-50:ATTGTCAGTCTGGTTAGAGCCACATGGCCTACCAGCTGACCAGTCCTTGCT

    // Ref ... GTGT                                                   TGTAG ...
    //             ATTGTCAGTCTGGTTAGAGCCACATGGCCTACCAGCTGACCAGTCCTTGCT
    //        
    // Alt ... GTGTATTGTCAGTCTGGTTAGAGCCACATGGCCTACCAGCTGACCAGTCCTTGCTTGTAG ...

    h.addVar(51, 50, "ATTGTCAGTCTGGTTAGAGCCACATGGCCTACCAGCTGACCAGTCCTTGCT");

    BOOST_CHECK_EQUAL(h.seq(0, 149), "TAATGACAGCGACTTGAGACATACACCTGCGCCTAATCACTTCAGAGGTGTATTGTCAGTCTGGTTAGAGCCACATGGCCTACCAGCTGACCAGTCCTTGCTTGTAGACGGGAAATAAGGATCTACCCTTACTTTGTGTCACTACTAGTGAAGCTGGCTCATGTGGGGACTAGGAAGAGCTGCCTCTCGTGGGAGCGATGC");

    int64_t min=-1, max=-1, sum=-1;
    h.lengthDiffs(min, max, sum);
    BOOST_CHECK_EQUAL(min, 51);
    BOOST_CHECK_EQUAL(max, 51);
    BOOST_CHECK_EQUAL(sum, 51);

}


BOOST_AUTO_TEST_CASE(haplotypeRange)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrQ", tp.string().c_str());

    h.addVar(5, 5, "CGCT");
    h.addVar(10, 10, "CA");

    // >chrQ
    // 0----5...----9.10---15---20
    // AAACCC---AAACC-CAAACCCGGGTTTGGGTTTGGGTTT
    //  [   |                   ]
    // AAACCCGCTAAACCAC...
    BOOST_CHECK_EQUAL(h.seq(2, 4), "ACC");
    BOOST_CHECK_EQUAL(h.seq(21, 25), "TTTGG");
    BOOST_CHECK_EQUAL(h.seq(2, 20), "ACCCGCTAAACCACAAACCCGGG");
    BOOST_CHECK_EQUAL(h.seq(6, 8), "GCT");
    // 8 is inside modified ref. We get everything from 8 in the modified region
    // 12 is outside, so we get everything up to the original 12
    // 
    BOOST_CHECK_EQUAL(h.seq(8, 12), "TAAACCACA");

    int64_t min=-1, max=-1, sum=-1;
    h.lengthDiffs(min, max, sum);
    BOOST_CHECK_EQUAL(min, 3);
    BOOST_CHECK_EQUAL(max, 4);
    BOOST_CHECK_EQUAL(sum, 4);
}

BOOST_AUTO_TEST_CASE(haplotypeRangeDel)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype h("chrQ", tp.string().c_str());

    h.addVar(5, 9, "AC");
    h.addVar(12, 13, "T");

    // >chrQ
    // AAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTT
    //      **     *
    // AAACCAC---CCT-ACCCGGGTTTGGGTTTGGGTTT
    BOOST_CHECK_EQUAL(h.seq(2, 4), "ACC");
    BOOST_CHECK_EQUAL(h.seq(21, 25), "TTTGG");
    BOOST_CHECK_EQUAL(h.seq(2, 15), "ACCACCCTAC");
    BOOST_CHECK_EQUAL(h.seq(6, 8), "CCC");
}
