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
 *  \brief Test Fasta file reading
 *
 * \file test_fasta.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/filesystem/path.hpp>

#include "Fasta.hh"

#include <iostream>

BOOST_AUTO_TEST_CASE(fastaRead)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    std::cerr << "Reading " << tp << std::endl;

    FastaFile f(tp.string().c_str());
    BOOST_CHECK_EQUAL(f.query("chrQ:5-9"), "CCAAA");
    BOOST_CHECK_EQUAL(f.query("chrT:10-33"), "TACAGGTCGTGTCGCGCAGAAGAA");
}


BOOST_AUTO_TEST_CASE(fastaReadEnd)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    std::cerr << "Reading " << tp << std::endl;

    FastaFile f(tp.string().c_str());
    BOOST_CHECK_EQUAL(f.query("chrS:151"), "");
}

BOOST_AUTO_TEST_CASE(fastaReadMultiline)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    std::cerr << "Reading " << tp << std::endl;

    FastaFile f(tp.string().c_str());
    BOOST_CHECK_EQUAL(f.query("chrS:20-147"), "CATACACCTGCGCCTAATCACTTCAGAGGTGTTGTAGACGGGAAATAAGGATCTACCCTTACTTTGTGTCACTACTAGTGAAGCTGGCTCATGTGGGGACTAGGAAGAGCTGCCTCTCGTGGGAGCGA");
}

BOOST_AUTO_TEST_CASE(fastaContigSizes)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                     .parent_path()   // test
                                     .parent_path()   // c++
                                 / boost::filesystem::path("data")
                                 / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.string().c_str());
    BOOST_CHECK_EQUAL(f.contigSize("chrT"), (size_t) 149);
    BOOST_CHECK_EQUAL(f.contigNonNSize("chrT"), (size_t) 149);
    BOOST_CHECK_EQUAL(f.contigSize("chrS"), (size_t) 150);
    BOOST_CHECK_EQUAL(f.contigNonNSize("chrS"), (size_t) 150);
    BOOST_CHECK_EQUAL(f.contigSize("chrU"), (size_t) 76);
    BOOST_CHECK_EQUAL(f.contigNonNSize("chrU"), (size_t) 36);
    BOOST_CHECK_EQUAL(f.contigSize("chrN"), (size_t) 76);
    BOOST_CHECK_EQUAL(f.contigNonNSize("chrN"), (size_t) 36);

    BOOST_CHECK_EQUAL(f.contigSize(), (size_t) 538);
    BOOST_CHECK_EQUAL(f.contigNonNSize(), (size_t) 458);
}
