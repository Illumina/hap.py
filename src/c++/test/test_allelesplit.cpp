// -*- mode: c++; indent-tabs-mode: nil; -*-
//
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
 *  \brief Test Allelic Primitive Splitting
 *
 *
 * \file test_allelesplit.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/filesystem/path.hpp>

#include "variant/VariantTee.hh"
#include "variant/VariantPrimitiveSplitter.hh"
#include "variant/VariantLocationAggregator.hh"
#include "variant/VariantAlleleUniq.hh"

#include <iostream>
#include <sstream>
#include <vector>
#include <set>

using namespace variant;

BOOST_AUTO_TEST_CASE(testVariantPrimitiveSplitter)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("chrQ.fa");

    VariantProcessor proc;
    VariantInjectorStep inj;
    proc.addStep(inj);

    VariantPrimitiveSplitter splitter;
    splitter.setReference(tp.c_str());
    proc.addStep(splitter);

    VariantLocationAggregator agg;
    agg.setAggregationType(variant::VariantLocationAggregator::AggregationType::aggregate_hetalt);
    proc.addStep(agg);

    VariantAlleleUniq uniq;
    proc.addStep(uniq);

    Variants vs;
    vs.chr = "chrS";
    vs.pos = 9;
    vs.len = 47 - 9 + 1;

    RefVar rv;
    rv.start = 9;
    rv.end = 40;
    // complex insertion
    rv.alt = "CGACTAGAGTCATACCCCCGCGCCTAATCACTTCACAGCTAATCACTAATCA";

    vs.variation.push_back(rv);

    rv.start = 9;
    rv.end = 47;
    // complex deletion
    rv.alt = "CGAACCGAGACATACAGCCTACTTCACAT";
    vs.variation.push_back(rv);

    vs.calls.resize(3);
    vs.calls[0].ngt = 2;
    vs.calls[0].gt[0] = 0;
    vs.calls[0].gt[1] = 1;

    vs.calls[1].ngt = 2;
    vs.calls[1].gt[0] = 0;
    vs.calls[1].gt[1] = 2;

    vs.calls[2].ngt = 2;
    vs.calls[2].gt[0] = 1;
    vs.calls[2].gt[1] = 2;

    inj.add(vs);

    // We have two complex alleles
    // ---------------------------
    //
    // (A)
    //
    // REF: CGACTTGAGACATACACCTGCGCCTAATCACT -------------------- TCAGAGG
    //           *   *     *  *              iiiiiiiiiiiiiiiiiiii
    // ALT: CGACTAGAGTCATACCCCCGCGCCTAATCACT TCACAGCTAATCACTAATCA TCAGAGG
    //
    //
    // (B)
    //
    // REF: CGACTTGAGACATACA CCTGC GCCTA ATCA CTTCAGAGG
    //         ***           ddddd       dddd      * *
    // ALT: CGAACCGAGACATACA ----- GCCTA ---- CTTCACAT-
    //
    //
    // Together, these alleles should come out:
    //
    // REF: CGACTTGAGACATACACCTGCGCCTAATCACT--------------------TCAGAGG
    //           *   *     *  *             iiiiiiiiiiiiiiiiiiii
    // ALTA CGACTAGAGTCATACCCCCGCGCCTAATCACTTCACAGCTAATCACTAATCATCAGAGG
    //         ***          ddddd     dddd                         * *d
    // ALTB CGAACCGAGACATACA-----GCCTA----CT--------------------TCACAT-
    //
    // We have GTs in three samples: SAMPLE1: 0/1  SAMPLE2: 0/2  SAMPLE3: 1/2

    std::map<int64_t, std::set<std::string> > expected_vars;

    // first three positions are homref, these are not reported

    // two het SNPs in 2 and 3
    expected_vars[12].insert("chrS:12-12 12-12:A . 0/1 0/1");
    expected_vars[13].insert("chrS:13-13 13-13:C . 0/1 0/1");
    expected_vars[14].insert("chrS:14-14 14-14:A 14-14:C 0/1 0/2 1/2");

    // homref GAG

    // A->T in 1 and 3
    expected_vars[18].insert("chrS:18-18 18-18:T 0/1 . 0/1");

    // homref CATAC

    // A->C in 1 and 3
    expected_vars[24].insert("chrS:24-24 24-24:C 0/1 . 0/1");
    // homref in 1 and hemi-homref in 3
    // het deletion in 2 and 3
    expected_vars[25].insert("chrS:25-29 25-29: . 0/1 0/1");

    // no-call in 2/3 because this is covered by het deletion

    // hemi-SNP in 3 (other allele is deletion) TODO: change this to report a 1 instead of 0/1
    expected_vars[27].insert("chrS:27-27 27-27:C 0/1 . 0/1");
    // homref GCCTA

    // deletion
    expected_vars[35].insert("chrS:35-38 35-38: . 0/1 0/1");

    // homref CT

    // het insertion in 1 and 3, homref in 2
    expected_vars[40].insert("chrS:40-40 41-40:TCACAGCTAATCACTAATCA 0/1 . 0/1");
    // homref CA

    // SNP homref SNP in 2 and 3
    expected_vars[44].insert("chrS:44-44 44-44:C . 0/1 0/1");
    expected_vars[46].insert("chrS:46-46 46-46:T . 0/1 0/1");

    // het deletion of length 1
    expected_vars[47].insert("chrS:47-47 47-47: . 0/1 0/1");

    std::set<int64_t> positions;
    while(proc.advance())
    {
        Variants & v = proc.current();
        std::cerr << v << "\n";
        positions.insert(v.pos);
        {
            std::ostringstream oss;
            oss << v;
            std::string res = stringutil::replaceAll(oss.str(), "1/0", "0/1");
            res = stringutil::replaceAll(res, "2/0", "0/2");
            res = stringutil::replaceAll(res, "2/1", "1/2");
            if(!expected_vars[v.pos].count(res))
            {
                std::cerr << "Unexpected variant calls: " << v << "\n";
                BOOST_ASSERT(false);
            }
        }
    }

    for(auto const & x : expected_vars)
    {
        if(!positions.count(x.first))
        {
            std::cerr << "Output for position " << x.first << " is missing: " << "\n";
            for(auto const & str : x.second)
            {
                std::cerr << "\t" << str << "\n";
            }
        }
    }
}
