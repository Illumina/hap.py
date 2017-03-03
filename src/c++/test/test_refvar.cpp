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
 *  \brief Test cases for haplotype comparison
 *
 *
 * \file test_refvar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>

#include "RefVar.hh"

// for testing realignRefVar
#include "Alignment.hh"

using namespace variant;

BOOST_AUTO_TEST_CASE(testRefVarLeftShiftSimple)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";

    leftShift(f, "chrQ", r);
    BOOST_CHECK_EQUAL(r.start, 17);
    BOOST_CHECK_EQUAL(r.end, 17);
    BOOST_CHECK_EQUAL(r.alt, "CGGGTTT");

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";

    // limit shift distance
    leftShift(f, "chrQ", r, 19);
    BOOST_CHECK_EQUAL(r.start, 19);
    BOOST_CHECK_EQUAL(r.end, 19);
    BOOST_CHECK_EQUAL(r.alt, "GGTTTGG");
}

BOOST_AUTO_TEST_CASE(testRefVarShiftNs)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 22;
    r.end = 23;
    r.alt = "C";

    leftShift(f, "chrU", r);
    BOOST_CHECK_EQUAL(r.start, 20);
    BOOST_CHECK_EQUAL(r.end, 21);
    BOOST_CHECK_EQUAL(r.alt, "A");

    r.start = 52;
    r.end = 52;
    r.alt = "GT";

    // limit shift distance
    rightShift(f, "chrU", r);
    BOOST_CHECK_EQUAL(r.start, 55);
    BOOST_CHECK_EQUAL(r.end, 55);
    BOOST_CHECK_EQUAL(r.alt, "TT");
}

BOOST_AUTO_TEST_CASE(testRefVarShift_HAP_64)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 23;
    r.end = 23;
    r.alt = "C";

    leftShift(f, "chrU", r);
    BOOST_CHECK_EQUAL(r.start, 23);
    BOOST_CHECK_EQUAL(r.end, 23);
    BOOST_CHECK_EQUAL(r.alt, "C");

    r.start = 23;
    r.end = 23;
    r.alt = "C";

    // limit shift distance
    rightShift(f, "chrU", r);
    BOOST_CHECK_EQUAL(r.start, 23);
    BOOST_CHECK_EQUAL(r.end, 23);
    BOOST_CHECK_EQUAL(r.alt, "C");
}

BOOST_AUTO_TEST_CASE(testRefVarLeftShiftList)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    std::list<RefVar> l;

    RefVar r;

    r.start = 5;
    r.end = 6;
    r.alt = "C";
    l.push_back(r);

    r.start = 8;
    r.end = 9;
    r.alt = "C";
    l.push_back(r);

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    l.push_back(r);

    leftShift<RefVar>(f, "chrQ", l);

    int count = 0;
    for(RefVar const & r : l)
    {
        switch(count++)
        {
            case 0:
                BOOST_CHECK_EQUAL(r.start, 5);
                BOOST_CHECK_EQUAL(r.end, 6);
                BOOST_CHECK_EQUAL(r.alt, "C");
                break;
            case 1:
                BOOST_CHECK_EQUAL(r.start, 7);
                BOOST_CHECK_EQUAL(r.end, 8);
                BOOST_CHECK_EQUAL(r.alt, "A");
                break;
            case 2:
                BOOST_CHECK_EQUAL(r.start, 17);
                BOOST_CHECK_EQUAL(r.end, 17);
                BOOST_CHECK_EQUAL(r.alt, "CGGGTTT");
                break;
            default:
                break;
        }
    }

}

BOOST_AUTO_TEST_CASE(testRefVarRightShiftSimple)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 5;
    r.end = 6;
    r.alt = "C";

    rightShift(f, "chrQ", r);
    BOOST_CHECK_EQUAL(r.start, 8);
    BOOST_CHECK_EQUAL(r.end, 9);
    BOOST_CHECK_EQUAL(r.alt, "C");

    // make sure we don't shift over the end
    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    rightShift(f, "chrQ", r);

    BOOST_CHECK_EQUAL(r.start, 35);
    BOOST_CHECK_EQUAL(r.end, 35);
    BOOST_CHECK_EQUAL(r.alt, "TGGGTTT");
}

BOOST_AUTO_TEST_CASE(testRefVarRightShiftList)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    std::list<RefVar> l;

    RefVar r;

    r.start = 5;
    r.end = 6;
    r.alt = "C";
    l.push_back(r);

    r.start = 7;
    r.end = 8;
    r.alt = "A";
    l.push_back(r);

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    l.push_back(r);

    rightShift<RefVar>(f, "chrQ", l);

    int count = 0;
    for(RefVar const & r : l)
    {
        switch(count++)
        {
            case 0:
                BOOST_CHECK_EQUAL(r.start, 6);
                BOOST_CHECK_EQUAL(r.end, 7);
                BOOST_CHECK_EQUAL(r.alt, "A");
                break;
            case 1:
                BOOST_CHECK_EQUAL(r.start, 8);
                BOOST_CHECK_EQUAL(r.end, 9);
                BOOST_CHECK_EQUAL(r.alt, "C");
                break;
            case 2:
                BOOST_CHECK_EQUAL(r.start, 35);
                BOOST_CHECK_EQUAL(r.end, 35);
                BOOST_CHECK_EQUAL(r.alt, "TGGGTTT");
                break;
            default:
                break;
        }
    }
}


BOOST_AUTO_TEST_CASE(testRefVarRightShiftMulti)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    // DEBUG *RV: 50-49:ATG
    // DEBUG *RV: 51-52:TG
    // DEBUG *RV: 56-57:GC
    // DEBUG *RV: 59-59:A
    // DEBUG *RV: 61-63:ACA
    // DEBUG *RV: 64-63:CCACA
    // DEBUG *RV: 67-67:C
    // DEBUG *RV: 69-70:
    // DEBUG *RV: 76-75:CAGCA
    // DEBUG *RV: 78-78:A
    // DEBUG *RV: 80-79:ACGAC
    // DEBUG *RV: 81-81:T
    // DEBUG *RV: 84-86:GTC

    RefVar r;
    std::list<RefVar> rl;
    r.start = 50;
    r.end = 49;
    r.alt = "ATG";
    rl.push_back(r);

    r.start = 51;
    r.end = 52;
    r.alt = "TG";
    rl.push_back(r);

    r.start = 56;
    r.end = 57;
    r.alt = "GC";
    rl.push_back(r);

    r.start = 59;
    r.end = 59;
    r.alt = "A";
    rl.push_back(r);

    r.start = 61;
    r.end = 63;
    r.alt = "ACA";
    rl.push_back(r);

    r.start = 64;
    r.end = 63;
    r.alt = "CCACA";
    rl.push_back(r);

    r.start = 67;
    r.end = 67;
    r.alt = "C";
    rl.push_back(r);

    r.start = 69;
    r.end = 70;
    r.alt = "";
    rl.push_back(r);

    r.start = 76;
    r.end = 75;
    r.alt = "CAGCA";
    rl.push_back(r);

    r.start = 78;
    r.end = 78;
    r.alt = "A";
    rl.push_back(r);

    r.start = 80;
    r.end = 79;
    r.alt = "ACGAC";
    rl.push_back(r);

    r.start = 81;
    r.end = 81;
    r.alt = "T";
    rl.push_back(r);

    r.start = 84;
    r.end = 86;
    r.alt = "GTC";
    rl.push_back(r);

    leftShift<RefVar>(f, "chrT", rl);

    const char * expected[14] = {
        "49-49:TATG",
        "51-52:TG",
        "56-57:GC",
        "59-59:A",
        "61-63:ACA",
        "64-63:CCACA",
        "67-67:C",
        "68-70:A",
        "75-75:CCAGCA",
        "78-78:A",
        "79-79:CACGAC",
        "81-81:T",
        "84-86:GTC"
    };

    int count = 0;
    for(RefVar const & rv : rl)
    {
        BOOST_CHECK(count < 14);
        std::ostringstream s;
        s << rv;
        BOOST_CHECK_EQUAL(s.str(), expected[count]);
        ++count;
    }
}

BOOST_AUTO_TEST_CASE(testRefVarAlleles)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    std::vector<RefVar> rl;

    // chrQ: AAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTT

    RefVar r;

    // long del
    r.start = 6;
    r.end = 11;
    r.alt = "A";
    rl.push_back(r);

    // short del
    r.start = 6;
    r.end = 9;
    r.alt = "A";
    rl.push_back(r);

    std::vector<std::string> alleles;
    toAlleles(f, "chrQ", rl, alleles);

    BOOST_CHECK_EQUAL(alleles.size(), (size_t)3);
    BOOST_CHECK_EQUAL("AAACCC", alleles[0]);
    BOOST_CHECK_EQUAL("A", alleles[1]);
    BOOST_CHECK_EQUAL("ACC", alleles[2]);
}

BOOST_AUTO_TEST_CASE(testAppendToVarList)
{
    std::list<RefVar> rvlist1;
    std::list<RefVar> rvlist2;
    std::list<RefVar> rvlist3;

    // Alignment matrix:
    //   0     .    :
    // REF AACACAC-AT
    //     ||||  | ||
    // PG  AATGTGCAAT
    //     ||||||  ||
    // QRY AACATGA-AT

    appendToVarList(1, 'A', 'A', rvlist1);
    appendToVarList(2, 'A', 'A', rvlist1);

    appendToVarList(3, 'C', 'T', rvlist2);
    appendToVarList(3, 'C', 'C', rvlist3);

    appendToVarList(4, 'A', 'G', rvlist2);
    appendToVarList(4, 'A', 'A', rvlist3);

    appendToVarList(5, 'C', 'T', rvlist1);
    appendToVarList(6, 'A', 'G', rvlist1);

    appendToVarList(7, 'C', 'C', rvlist2);
    appendToVarList(7, 'C', 'A', rvlist3);

    appendToVarList(8, '-', 'A', rvlist2);
    appendToVarList(8, '-', '-', rvlist3);

    appendToVarList(9,  'A', 'A', rvlist1);
    appendToVarList(10, 'T', 'T', rvlist1);

    // for(auto & x : rvlist1)
    // {
    //     std::cout << x << "; ";
    // }
    // std::cout << "\n";

    // for(auto & x : rvlist2)
    // {
    //     std::cout << x << "; ";
    // }
    // std::cout << "\n";

    // for(auto & x : rvlist3)
    // {
    //     std::cout << x << "; ";
    // }
    // std::cout << "\n";

    BOOST_CHECK_EQUAL(rvlist1.size(), (size_t)1);
    BOOST_CHECK_EQUAL(rvlist2.size(), (size_t)2);
    BOOST_CHECK_EQUAL(rvlist3.size(), (size_t)1);

    BOOST_CHECK_EQUAL(rvlist1.front().repr(), "5-6:TG");
    BOOST_CHECK_EQUAL(rvlist2.front().repr(), "3-4:TG");
    BOOST_CHECK_EQUAL(rvlist2.back().repr(), "8-7:A");
    BOOST_CHECK_EQUAL(rvlist3.front().repr(), "7-7:A");
}

BOOST_AUTO_TEST_CASE(testAppendToVarList2)
{
    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG
    //                                                    *        *  *                                                                       *   *
    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG------------------------------------------------------------------TTTTGC---------G

    std::string s1 = "GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG";
    std::string s2 = "GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG------------------------------------------------------------------TTTTGC---------G";

    std::list<RefVar> rvlist;

    for (size_t i = 0; i < s1.size(); ++i)
    {
        appendToVarList(i+1,  s1[i], s2[i], rvlist);
    }

    std::ostringstream oss;
    for(auto & x : rvlist)
    {
        oss << x << "; ";
    }
    BOOST_CHECK_EQUAL(oss.str(), "52-52:A; 54-59:; 61-61:C; 64-64:C; 69-134:; 136-136:T; 140-149:C; ");
}

BOOST_AUTO_TEST_CASE(testAppendToVarList3)
{
    std::string s1 = "GAAGTACAGAGTCGATTTG--GACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCTTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG";
    std::string s2 = "GAAGTACAGAGTCGATTTGGA--CGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG------------------------------------------------------------------TTTTGC---------G";

    std::list<RefVar> rvlist;

    size_t refpos = 0;
    for (size_t i = 0; i < s1.size(); ++i)
    {
        appendToVarList(refpos+1, s1[i], s2[i], rvlist);
        if(s1[i] != '-')
        {
            ++refpos;
        }
    }

    std::ostringstream oss;
    for(auto & x : rvlist)
    {
        oss << x << "; ";
    }
    BOOST_CHECK_EQUAL(oss.str(), "20-21:GA; 52-52:A; 54-61:CC; 64-64:C; 69-134:; 136-136:T; 140-149:C; ");
}

BOOST_AUTO_TEST_CASE(testAppendToVarList4)
{
    std::string s1 = "TG--T";
    std::string s2 = "TGGAT";

    std::list<RefVar> rvlist;

    size_t refpos = 0;
    for (size_t i = 0; i < s1.size(); ++i)
    {
        appendToVarList(refpos+1, s1[i], s2[i], rvlist);
        if(s1[i] != '-')
        {
            ++refpos;
        }
    }

    std::ostringstream oss;
    for(auto & x : rvlist)
    {
        oss << x << "; ";
    }
    BOOST_CHECK_EQUAL(oss.str(), "3-2:GA; ");
}

BOOST_AUTO_TEST_CASE(testRefVarApply)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    RefVar rv;
    int64_t s, e;

    rv.start = 8;
    rv.end = 7;
    rv.alt = "AAACCC";
    s =  0;
    e = 15;

    std::string res = rv.apply(f, "chrT", s, e);

    // ref: ATTCTGAC------ATACAGGT
    // alt: ATTCTGACaaacccATACAGGT

    BOOST_CHECK_EQUAL(res, "ATTCTGACAAACCCATACAGGT");

    rv.start = 8;
    rv.end = 10;
    rv.alt = "AAACCC";
    s =  0;
    e = 15;

    // ref: ATTCTGAC------ACAGGT
    // alt: ATTCTGACaaacccACAGGT

    res = rv.apply(f, "chrT", s, e);
    BOOST_CHECK_EQUAL(res, "ATTCTGACAAACCCCAGGT");
}


BOOST_AUTO_TEST_CASE(testRefVarPrimitives)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());

    //>chrS:10-15
    //CGACT

    RefVar rv;
    rv.start = 9;
    rv.end = 14;
    // complex insertion
    rv.alt = "TGCCTTT";

    std::list<RefVar> rvl;
    toPrimitives(f, "chrS", rv, rvl);

    {
        std::ostringstream oss;
        for(auto & x : rvl)
        {
            oss << x << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "9-9:T; 11-11:C; 15-14:T; ");
    }

    rv.start = 9;
    rv.end = 14;
    // complex deletion
    rv.alt = "TGCC";

    rvl.clear();
    toPrimitives(f, "chrS", rv, rvl);

    {
        std::ostringstream oss;
        for(auto & x : rvl)
        {
            oss << x << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "9-9:T; 11-11:C; 13-14:; ");
    }
}

BOOST_AUTO_TEST_CASE(testRefVarPrimitiveAlign)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());
    Alignment * aln = makeAlignment("klibg");

    //>chrS
    //CGACT

    RefVar rv;
    rv.start = 9;
    rv.end = 14;
    // complex insertion
    rv.alt = "TGCCTTT";

    std::list<RefVar> rvl;
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for(auto & x : rvl)
        {
            oss << x << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "9-9:T; 11-11:C; 15-14:T; ");
    }

    rv.start = 9;
    rv.end = 14;
    // complex deletion
    rv.alt = "TGCC";

    rvl.clear();
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for(auto & x : rvl)
        {
            oss << x << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "9-9:T; 11-11:C; 13-14:; ");
    }
    delete aln;
}

BOOST_AUTO_TEST_CASE(testRefVarPrimitiveAlign2)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    FastaFile f(tp.c_str());
    Alignment * aln = makeAlignment("klibg");

    //>chrS:10
    //CGACTTGAGACATACACCTGCGCCTAATCACTTCAGAG...

    RefVar rv;
    rv.start = 9;
    rv.end = 40;
    // complex insertion
    rv.alt = "CGACTAGAGTCATACCCCCGCGCCTAATCACTTCACAGCTAATCACTAATCA";

    std::list<RefVar> rvl;
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for(auto & x : rvl)
        {
            oss << x << "; ";
        }
        // REF: CGACTTGAGACATACACCTGCGCCTAATCACTTCAGAG--------------
        //           *   *     *  *                   iiiiiiiiiiiiii
        // ALT: CGACTAGAGTCATACCCCCGCGCCTAATCACTTCACAGCTAATCACTAATCA
        BOOST_CHECK_EQUAL(oss.str(), "14-14:A; 18-18:T; 24-24:C; 27-27:C; 41-40:TCACAGCTAATCACTAATCA; ");
    }

    rv.start = 9;
    rv.end = 47;
    // complex deletion
    rv.alt = "CGAACCGAGACATACAGCCTACTTCACAT";

    rvl.clear();
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for(auto & x : rvl)
        {
            oss << x << "; ";
        }
        // REF: CGACTTGAGACATACA CCTGC GCCTA ATCA CTTCAGAGG
        //         ***           ddddd       dddd      * *
        // ALT: CGAACCGAGACATACA ----- GCCTA ---- CTTCACAT-
        BOOST_CHECK_EQUAL(oss.str(), "12-12:A; 13-13:C; 14-14:C; 25-29:; 35-38:; 44-44:C; 46-46:T; 47-47:; ");
    }
    delete aln;
}

BOOST_AUTO_TEST_CASE(testRefVarLeftShiftStartOfChr)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                         .parent_path()   // test
                                         .parent_path()   // c++
                                 / boost::filesystem::path("data")
                                 / boost::filesystem::path("leftshifting_example")
                                 / boost::filesystem::path("ref.fa");

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 2;
    r.end = 3;
    r.alt = "";

    leftShift(f, "chrT", r);
    BOOST_CHECK_EQUAL(r.start, 0);
    BOOST_CHECK_EQUAL(r.end, 2);
    BOOST_CHECK_EQUAL(r.alt, "C");

    r.start = 2;
    r.end = 3;
    r.alt = "";

    leftShift(f, "chrT", r, -1000, true);
    BOOST_CHECK_EQUAL(r.start, 0);
    BOOST_CHECK_EQUAL(r.end, 2);
    BOOST_CHECK_EQUAL(r.alt, "C");
}

