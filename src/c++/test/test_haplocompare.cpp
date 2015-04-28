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
 * \file test_haplocompare.cpp
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

#include "Haplotype.hh"
#include "HaploCompare.hh"

using namespace variant;
using namespace haplotypes;
 
BOOST_AUTO_TEST_CASE(haplocompareBasic)
{
    HaploCompare hc;

    // CTGTGTGTG------TGTGTGAAAA
    // CTGTGTGTGAGTGTGTGTGTGAAAA

    hc.setRef("CTGTGTGTGTGTGTGAAAA");
    hc.setAlt("CTGTGTGTGAGTGTGTGTGTGAAAA");

    std::list<RefVar> vlist;
    hc.getVariants(vlist);
    Alignment * al = hc.getAlignment();

    int r0, r1, a0, a1;
    std::string cig;
    al->getCigar(r0, r1, a0, a1, cig);

    BOOST_CHECK_EQUAL(r0, 0);
    BOOST_CHECK_EQUAL(r1, 18);
    BOOST_CHECK_EQUAL(a0, 0);
    BOOST_CHECK_EQUAL(a1, 24);
    BOOST_CHECK_EQUAL(cig, "9M6D10M");

    BOOST_CHECK_EQUAL(vlist.size(), (size_t)1);

    BOOST_CHECK_EQUAL(vlist.front().start, 9);
    BOOST_CHECK_EQUAL(vlist.front().end, 8);
    BOOST_CHECK_EQUAL(vlist.front().alt, "AGTGTG");

    // CTGTGTGTGTGTGTGTGTGTGAAAA
    // CTGTGTGTGTGTGT------GAAAA
    hc.setRef("CTGTGTGTGTGTGTGTGTGTGAAAA");
    hc.setAlt("CTGTGTGTGTGTGTGAAAA");

    vlist.clear();
    hc.getVariants(vlist);
    al->getCigar(r0, r1, a0, a1, cig);

    BOOST_CHECK_EQUAL(r0, 0);
    BOOST_CHECK_EQUAL(r1, 24);
    BOOST_CHECK_EQUAL(a0, 0);
    BOOST_CHECK_EQUAL(a1, 18);
    BOOST_CHECK_EQUAL(cig, "15M6I4M");

    BOOST_CHECK_EQUAL(vlist.size(), (size_t)1);

    BOOST_CHECK_EQUAL(vlist.front().start, 15);
    BOOST_CHECK_EQUAL(vlist.front().end, 20);
    BOOST_CHECK_EQUAL(vlist.front().alt, "");

    // CTGTGTGTGTGTGTGTGTGTGAAAA
    //          *
    // CTGTGTGTGAGTGTGTGTGTGAAAA
    hc.setRef("CTGTGTGTGTGTGTGTGTGTGAAAA");
    hc.setAlt("CTGTGTGTGAGTGTGTGTGTGAAAA");

    vlist.clear();
    hc.getVariants(vlist);
    al->getCigar(r0, r1, a0, a1, cig);

    BOOST_CHECK_EQUAL(r0, 0);
    BOOST_CHECK_EQUAL(r1, 24);
    BOOST_CHECK_EQUAL(a0, 0);
    BOOST_CHECK_EQUAL(a1, 24);
    BOOST_CHECK_EQUAL(cig, "25M");

    BOOST_CHECK_EQUAL(vlist.size(), (size_t)1);

    BOOST_CHECK_EQUAL(vlist.front().start, 9);
    BOOST_CHECK_EQUAL(vlist.front().end, 9);
    BOOST_CHECK_EQUAL(vlist.front().alt, "A");

}

BOOST_AUTO_TEST_CASE(haplocompareComplex)
{
    HaploCompare hc;

    // platinum complex variant at chr1:17678525
    // ATGGGAGGAGAGGCGGCCAGGGAGGTGATGGGAGGAGAGGCAGTCAGGGAG
    //                          TGATGGG                          
    //                          AGACAGGAGGCAGT

    hc.setRef("ATGGGAGGAGAGGCGGCCAGGGAGGTGATGGGAGGAGAGGCAGTCAGGGAG");
    hc.setAlt("ATGGGAGGAGAGGCGGCCAGGGAGGAGACAGGAGGCAGTAGGAGAGGCAGTCAGGGAG");

    std::list<RefVar> vlist;
    hc.getVariants(vlist);
    Alignment * al = hc.getAlignment();

    int r0, r1, a0, a1;
    std::string cig;
    al->getCigar(r0, r1, a0, a1, cig);

    BOOST_CHECK_EQUAL(r0, 0);
    BOOST_CHECK_EQUAL(r1, 50);
    BOOST_CHECK_EQUAL(a0, 0);
    BOOST_CHECK_EQUAL(a1, 57);
    BOOST_CHECK_EQUAL(cig, "35M7D16M");

    // ATGGGAGGAGAGGCGGCCAGGGAGGTGATGGGAGG-------AGAGGCAGTCAGGGAG
    //                          *  **            
    // ATGGGAGGAGAGGCGGCCAGGGAGGAGACAGGAGGCAGTAGGAGAGGCAGTCAGGGAG

    BOOST_CHECK_EQUAL(vlist.size(), (size_t)3);

    int count = 0;
    for (RefVar const & rv : vlist)
    {
        switch(count++)
        {
            case 0:
                BOOST_CHECK_EQUAL(rv.start, 25);
                BOOST_CHECK_EQUAL(rv.end, 25);
                BOOST_CHECK_EQUAL(rv.alt, "A");
                break;
            case 1:
                BOOST_CHECK_EQUAL(rv.start, 28);
                BOOST_CHECK_EQUAL(rv.end, 29);
                BOOST_CHECK_EQUAL(rv.alt, "CA");
                break;
            case 2:
                BOOST_CHECK_EQUAL(rv.start, 35);
                BOOST_CHECK_EQUAL(rv.end, 34);
                BOOST_CHECK_EQUAL(rv.alt, "CAGTAGG");
                break;
            default:
                BOOST_CHECK(false);
        }
    }

    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data") 
                                    / boost::filesystem::path("chrQ.fa");

    Haplotype ht("chrR", tp.string().c_str());
    ht.addVar(vlist);
    BOOST_CHECK_EQUAL(ht.seq(0, 50), "ATGGGAGGAGAGGCGGCCAGGGAGGAGACAGGAGGCAGTAGGAGAGGCAGTCAGGGAG");
}

BOOST_AUTO_TEST_CASE(haplocompareRandom)
{
    const int R = 150, M = 50, 
#ifdef _DEBUG
    N = 50;
#else
    N = 1000;
#endif

    for (int j = 0; j < N; ++j) 
    {
        HaploCompare hc;
        // make a reference sequence
        boost::filesystem::path temp = boost::filesystem::unique_path("%%%%-%%%%-%%%%-%%%%.fa");
        std::ostringstream sref;

        const char * nts[4] = {"A", "C", "G", "T"};
        for (int i = 0; i < R; ++i)
        {
            const char * r = nts[rand() & 3];
            sref << r;
        }
        std::string ref(sref.str());

        std::ofstream rseq(temp.native());
        rseq << ">chrQ" << "\n";
        rseq << ref << "\n";
        rseq.close();

        std::string alt(ref);
        // make an alt sequence
        for (int i = M; (i < 2*M && i < (int)alt.size()); ++i)
        {
            switch(rand() % 3)
            {
                case 0: // randomly change
                    alt[i] = nts[rand() & 3][0];
                    break;
                case 1: // random insertion
                {
                    std::string random;
                    size_t rlen = rand() % 20;
                    for(size_t j = 0; j < rlen; ++j)
                    {
                        random += nts[rand() & 3][0];
                    }
                    alt.insert(i, random);
                    i += rlen-1;
                    break;
                }
                case 2: // random deletion
                {
                    std::string random;
                    size_t rlen = rand() % 20;
                    if(i + rlen > alt.size())
                    {
                        rlen = alt.size() - i;
                    }
                    alt.replace(i, rlen, "");
                    break;
                }
            }
        }

        hc.setRef(ref.c_str());
        hc.setAlt(alt.c_str());
        std::list<RefVar> vlist;
        hc.getVariants(vlist);

        {
            Haplotype ht("chrQ", temp.c_str());
            std::string halt = "-";
            try
            {
                ht.addVar(vlist);
                halt = ht.seq(0, ref.size()-1);
            }
            catch(std::runtime_error & e)
            {
                std::cerr << "Error creating the haplotype sequence: " << e.what() << "\n";
            }

            BOOST_CHECK_EQUAL(alt, halt);

            if(alt != halt)
            {
                Alignment * aln = hc.getAlignment();
                int r0, r1, a0, a1;
                std::string cig;
                aln->getCigar(r0, r1, a0, a1, cig);
                std::cerr << "DEBUG REF: " << ref << "\n";
                std::cerr << "DEBUG ALT: " << alt << "\n";
                std::cerr << "DEBUG ALN: " << r0 << " " << r1 << " " << a0 << " " << a1 << " : " << cig << "\n";
                for (RefVar const & rv : vlist)
                {
                    std::cerr << "DEBUG RV: " << rv << "\n";
                }
                break;
            }            
        }

        std::list<RefVar> vlist_b;
        // same with left-shifting afterwards
        {
            Haplotype ht("chrQ", temp.c_str());
            FastaFile f(temp.c_str());
            vlist_b = vlist;
            leftShift(f, "chrQ", vlist);

            std::string halt = "-";
            try
            {
                ht.addVar(vlist);
                halt = ht.seq(0, ref.size()-1);
            }
            catch(std::runtime_error & e)
            {
                std::cerr << "Error creating the haplotype sequence: " << e.what() << "\n";
            }

            BOOST_CHECK_EQUAL(alt, halt);

            if(alt != halt)
            {
                Alignment * aln = hc.getAlignment();
                int r0, r1, a0, a1;
                std::string cig;
                aln->getCigar(r0, r1, a0, a1, cig);
                std::cerr << "DEBUG REF: " << ref << "\n";
                std::cerr << "DEBUG ALT: " << alt << "\n";
                std::cerr << "DEBUG ALN: " << r0 << " " << r1 << " " << a0 << " " << a1 << " : " << cig << "\n";
                for (RefVar const & rv : vlist_b)
                {
                    std::cerr << "DEBUG *RV: " << rv << "\n";
                }

                for (RefVar const & rv : vlist)
                {
                    std::cerr << "DEBUG <RV: " << rv << "\n";
                }
                break;
            }            
        }

        // same with right-shifting afterwards
        {
            Haplotype ht("chrQ", temp.c_str());
            FastaFile f(temp.c_str());
            vlist = vlist_b;
            rightShift<RefVar>(f, "chrQ", vlist);

            std::string halt = "-";
            try
            {
                ht.addVar(vlist);
                halt = ht.seq(0, ref.size()-1);
            }
            catch(std::runtime_error & e)
            {
                std::cerr << "Error creating the haplotype sequence: " << e.what() << "\n";
            }

            BOOST_CHECK_EQUAL(alt, halt);

            if(alt != halt)
            {
                Alignment * aln = hc.getAlignment();
                int r0, r1, a0, a1;
                std::string cig;
                aln->getCigar(r0, r1, a0, a1, cig);
                std::cerr << "DEBUG REF: " << ref << "\n";
                std::cerr << "DEBUG ALT: " << alt << "\n";
                std::cerr << "DEBUG ALN: " << r0 << " " << r1 << " " << a0 << " " << a1 << " : " << cig << "\n";

                for (RefVar const & rv : vlist_b)
                {
                    std::cerr << "DEBUG *RV: " << rv << "\n";
                }

                for (RefVar const & rv : vlist)
                {
                    std::cerr << "DEBUG >RV: " << rv << "\n";
                }
                break;
            }            
        }

        std::cerr << ".";
        Haplotype::resetRefs();

        boost::filesystem::remove(temp);
        temp += ".fai";
        boost::filesystem::remove(temp);
    }
    std::cerr << "\n";
}
