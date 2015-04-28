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
 * Test for SharedVar class
 * 
 * \file test_sharedvar.cpp
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

#include "SharedVar.hh"
#include "helpers/Timing.hh"

#include "Error.hh"

// #define DEBUG_SHAREDVAR

using namespace variant;
using namespace haplotypes;
 
#ifdef HAS_MUSCLE

BOOST_AUTO_TEST_CASE(sharedVarBasic)
{
    SharedVar sv;

    std::list<RefVar> vars;

    sv.splitVar(
        0,
        "AACACACAT",    
        "AACATGCAAT",   
        "AACATGAAT",    
        vars
    );

    // AACACAC-AT
    // AACATGCAAT
    // AACATGA-AT

    {
        std::ostringstream oss;
        for(auto & x : vars)
        {
            oss << x << ":" << x.flags << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "4-5:TG:3; 6-6:A:2; 7-6:A:1; ");
    }
}


BOOST_AUTO_TEST_CASE(sharedVar2)
{
    SharedVar sv;

    std::list<RefVar> vars;


    // TGACAGGATACATACATACATACATACATACATAC--------------ATATATATATATATATATATATATATAAAG
    // TGACAGGATACATACATACATACATACATACATACATACATATATATATATATATATATATATATATATATATATAAAG
    // TGACAGGATACATACATACATACATACATACATAT--------ACATACATATATATATATATATATATATATATAAAG

    sv.splitVar(
        0,
        "TGACAGGATACATACATACATACATACATACATACATATATATATATATATATATATATATAAAG",    
        "TGACAGGATACATACATACATACATACATACATACATACATATATATATATATATATATATATATATATATATATAAAG",   
        "TGACAGGATACATACATACATACATACATACATATACATACATATATATATATATATATATATATATAAAG",
        vars
    );

    {
        std::ostringstream oss;
        for(auto & x : vars)
        {
            oss << x << ":" << x.flags << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "34-34:T:2; 35-34:ATACATAT:1; 35-34:A:3; 35-34:T:1; 35-34:C:2; 35-34:ATA:3; 35-34:T:1; 35-34:C:2; ");
    }
}

BOOST_AUTO_TEST_CASE(sharedVar3)
{
    SharedVar sv;

    std::list<RefVar> vars;

    sv.splitVar(
        0,
        "GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG",    
        "GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG",   
        "GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG",
        vars
    );

    // GGTAATTCTAAGCAAAAACACACACACCAGA----ATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG
    // GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG
    // GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG

    {
        std::ostringstream oss;
        for(auto & x : vars)
        {
            oss << x << ":" << x.flags << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "31-30:ATAT:3; ");
    }
}

BOOST_AUTO_TEST_CASE(sharedVar4)
{
    SharedVar sv;

    std::list<RefVar> vars;

    sv.splitVar(
        0,
        "GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG",    
        "GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG",   
        "GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG",
        vars
    );

    // GGTAATTCTAAGCAAAAACACACACACCAGA------ATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG
    // GGTAATTCTAAGCAAAAACACACACACCAGA------ATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG
    // GGTAATTCTAAGCAAAAACACACACACCAGAATATATATATATATATATATATATATATAAATCACTGTAAATACATATAAACTG

    {
        std::ostringstream oss;
        for(auto & x : vars)
        {
            oss << x << ":" << x.flags << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "31-30:ATATAT:2; ");
    }
}

BOOST_AUTO_TEST_CASE(sharedVar5)
{
    SharedVar sv;

    std::list<RefVar> vars;

    sv.splitVar(
        0,
        "GACTGAGGGAGACTCCGTCTCAAGAAAAAAAAAAAAAAAAAAAAAAAGAAATGT",    
        "GACAGAGCAAGACTCCGTCTCAAGAAAAAAAAAAAAAAAAAAAAAAAGAAATGT",   
        "GACTGAGGGAGACTCCGTCTCAAAAATCAAGAAAAAAAAAAAAAAAAAAAAAAAGAAATGT",
        vars
    );

    // GACTGAGGGAGACTCCGTC-------TCAAGAAAAAAAAAAAAAAAAAAAAAAAGAAATGT
    // GACAGAGCAAGACTCCGTC-------TCAAGAAAAAAAAAAAAAAAAAAAAAAAGAAATGT
    // GACTGAGGGAGACTCCGTCTCAAAAATCAAGAAAAAAAAAAAAAAAAAAAAAAAGAAATGT

    {
        std::ostringstream oss;
        for(auto & x : vars)
        {
            oss << x << ":" << x.flags << "; ";
        }
        BOOST_CHECK_EQUAL(oss.str(), "3-3:A:1; 7-8:CA:1; 19-18:TCAAAAA:2; ");
    }
}

struct SharedVarTimer
{
    SharedVarTimer(size_t len)
    {
        for(size_t i = 0; i < len; ++i)
        {
            ref += chars[rand() & 3];

            int rp = rand() & 10;

            if(rp < 5)
            {
                alt1 += ref.back();
                alt2 += ref.back();
            }
            else if(rp < 8)
            {
                alt1 += chars[rand() & 3];
                alt2 += chars[rand() & 3];
            }
            else
            {
                rp = rand() & 4;
                if (rp < 1)
                {
                    alt1 += chars[rand() & 3];
                }
                else if (rp < 2)
                {
                    alt2 += chars[rand() & 3];
                }
                else
                {
                    ref += chars[rand() & 3];
                }
            }
        }
    }

    void operator() ()
    {
        std::list<RefVar> vars;

        sv.splitVar(
            0,
            ref,    
            alt1,   
            alt2,
            vars
        );
    }

    SharedVar sv;
    std::string ref;
    std::string alt1, alt2;

    const char chars[4] = {'A', 'C', 'G', 'T'};
};

BOOST_AUTO_TEST_CASE(sharedVarPerformance)
{
    std::cerr << "Testing Shared Variant Multi-aligner performance." << "\n";
    std::cerr << "\n";
    int count = 0;

    const int SIZE = 100;
#ifdef _DEBUG
    const int N = 1;
    const int REPS = 2;
#else
    const int N = 100;    
    const int REPS = 32;
#endif
    while(count++ < REPS)
    {
        double result;
        TIMEIT(result, N, SharedVarTimer, SIZE);
        std::cerr << result << "\n";
    }
    std::cerr << "\n";

}

struct SharedVarTester
{
    SharedVarTester(size_t len)
    {
        for(size_t i = 0; i < len; ++i)
        {
            ref += chars[rand() & 3];

            int rp = rand() & 10;

            if(rp < 5)
            {
                alt1 += ref.back();
                alt2 += ref.back();
            }
            else if(rp < 8)
            {
                alt1 += chars[rand() & 3];
                alt2 += chars[rand() & 3];
            }
            else
            {
                rp = rand() & 4;
                if (rp < 1)
                {
                    alt1 += chars[rand() & 3];
                }
                else if (rp < 2)
                {
                    alt2 += chars[rand() & 3];
                }
                else
                {
                    ref += chars[rand() & 3];
                }
            }
        }
    }

    void operator() ()
    {
        std::list<RefVar> vars;

        sv.splitVar(
            0,
            ref,    
            alt1,   
            alt2,
            vars
        );

        auto applyVar = [](std::string & s, RefVar & rv, int64_t & shift)
        {
            int64_t pos = rv.start - shift;
            int64_t reflen = rv.end - rv.start + 1;
            int64_t altlen = (int64_t)rv.alt.size();

            assert(reflen >= 0);

            if(reflen == 0)
            {
                if(pos < (signed)s.size())
                {
                    s.insert(pos, rv.alt);
                }
                else if(pos == (signed)s.size())
                {
                    s+= rv.alt;
                }
                else
                {
                    error("Pos %i out of range for s=%s (l=%i)", pos, s.c_str(), s.size());
                }
            }
            else
            {
                s.replace(pos, reflen, rv.alt);
            }

            shift += reflen - altlen;
        };


        std::string aref1 = ref;
        auto sp = vars.begin();
        int64_t s_ap1 = 0;

        while(sp != vars.end())
        {
            if (sp->flags == 1 || sp->flags == 3)
            {
#ifdef DEBUG_SHAREDVAR
                std::cerr << "s1:" << aref1 << "[" << s_ap1 << "]" << " + " << *sp << "\n";
#endif                
                applyVar(aref1, *sp, s_ap1);
            }
            ++sp;
        }

        BOOST_CHECK_EQUAL(alt1, aref1);

        std::string aref2 = ref;
        sp = vars.begin();
        int64_t s_ap2 = 0;

        while(sp != vars.end())
        {
            // process by end pos
            if (sp->flags == 2 || sp->flags == 3)
            {
                // support insertions with no reference character
#ifdef DEBUG_SHAREDVAR
                std::cerr << "s2:" << aref2 << "[" << s_ap2 << "]" << " + " << *sp << "\n";
#endif                
                applyVar(aref2, *sp, s_ap2);
            }
            ++sp;
        }
        BOOST_CHECK_EQUAL(alt2, aref2);
    }

    SharedVar sv;
    std::string ref;
    std::string alt1, alt2;

    const char chars[4] = {'A', 'C', 'G', 'T'};
};

BOOST_AUTO_TEST_CASE(sharedVarConsistency)
{
    std::cerr << "Testing Shared Variant Multi-aligner." << "\n";
    std::cerr << "\n";
    int count = 0;

    auto nop = [](double) {};

    const int SIZE = 100;
#ifdef _DEBUG
    const int N = 10;
    const int REPS = 1;
#else
    const int N = 100;    
    const int REPS = 32;
#endif
    while(count++ < REPS)
    {
        double result;
        TIMEIT(result, N, SharedVarTester, SIZE);
        // avoids a warning
        nop(result);
    }
    std::cerr << "\n";

}

#endif

