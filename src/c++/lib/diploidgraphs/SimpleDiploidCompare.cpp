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
 * Simple Diploid comparison by variant position and count
 *
 * \file SimpleDiploidCompare.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "SimpleDiploidCompare.hh"

// #define DEBUG_SIMPLECMP

#include "Variant.hh"
#include "RefVar.hh"
#include "Haplotype.hh"

#include <map>
#include <set>
#include <cmath>
#include <cassert>
#ifdef DEBUG_SIMPLECMP
#include <fstream>
#endif

namespace haplotypes
{

/**
 * @brief Compare two variants
 *
 */
DiploidComparisonOutcome compareVariants(variant::Variants & vars, int r1, int r2, int & non_snp_alleles, int & calls_1, int & calls_2,
                                         bool use_filtered_1, bool use_filtered_2)
{
    using namespace variant;
#ifdef DEBUG_SIMPLECMP
    std::cerr << "SimpleDiploidCompare: " << vars << "\n";
#endif
    Call & call1 = vars.calls[r1];
    Call & call2 = vars.calls[r2];
    gttype gtt1 = getGTType(call1);
    gttype gtt2 = getGTType(call2);

#ifdef DEBUG_SIMPLECMP
    std::cerr << "                Call1:" << call1 << " (" << gtt1 << ")" << "\n";
    std::cerr << "                Call2:" << call2 << " (" << gtt2 << ")" << "\n";
#endif

    bool is_filtered_1 = false;
    if(!use_filtered_1)
    {
        for(size_t f = 0; f < call1.nfilter; ++f)
        {
            if(call1.filter[f] != "." && call1.filter[f] != "PASS")
            {
                is_filtered_1 = true;
                break;
            }
        }
    }

    bool is_filtered_2 = false;
    if(!use_filtered_2)
    {
        for(size_t f = 0; f < call2.nfilter; ++f)
        {
            if(call2.filter[f] != "." && call2.filter[f] != "PASS")
            {
                is_filtered_2 = true;
                break;
            }
        }
    }

    bool is_match = false;
    const bool is_ambig_1 = int(vars.ambiguous_alleles.size()) > r1 && !vars.ambiguous_alleles[r1].empty();
    const bool is_ambig_2 = int(vars.ambiguous_alleles.size()) > r2 && !vars.ambiguous_alleles[r2].empty();

    const bool is_call_1 = !(gtt1 == gt_homref || gtt1 == gt_unknown || is_ambig_1 || is_filtered_1);
    const bool is_call_2 = !(gtt2 == gt_homref || gtt2 == gt_unknown || is_ambig_2 || is_filtered_2);

    if (is_call_1)
    {
        ++calls_1;
    }

    if (is_call_2)
    {
        ++calls_2;
    }

    if(is_ambig_1 || is_ambig_2)
    {
        vars.setInfo("type", "N");
        vars.setInfo("kind", "ambiguous");
        return dco_mismatch;
    }

    static const char * gttype_string[] = {
        "gt_homref",
        "gt_haploid",
        "gt_het",
        "gt_homalt",
        "gt_hetalt",
        "gt_unknown"
    };

    if (is_call_1)
    {
        size_t ngt = call1.ngt;
        if(gtt1 == gt_homalt)
        {
            // only write one record for homalts
            ngt = 1;
        }
        for (size_t i = 0; i < ngt; ++i)
        {
            if(call1.gt[i] < 1)
            {
                continue;
            }
            RefVar & rv (vars.variation[call1.gt[i]-1]);
            int64_t reflen = rv.end - rv.start + 1;
            int64_t altlen = (int64_t)(rv.alt.size());

            if (reflen != 1 || altlen != 1)
            {
                ++non_snp_alleles;
            }
        }
    }

    if (is_call_2)
    {
        size_t ngt = call2.ngt;
        if(gtt2 == gt_homalt)
        {
            // only write one record for homalts
            ngt = 1;
        }
        for (size_t i = 0; i < ngt; ++i)
        {
            if(call2.gt[i] < 1)
            {
                continue;
            }
            RefVar & rv (vars.variation[call2.gt[i]-1]);
            int64_t reflen = rv.end - rv.start + 1;
            int64_t altlen = (int64_t)(rv.alt.size());

            if (reflen != 1 || altlen != 1)
            {
                ++non_snp_alleles;
            }
        }
    }

    // only 1?
    if (is_call_1 && !is_call_2)
    {
        // missing from "1" = FN
        vars.setInfo("type", "FN");
        vars.setInfo("kind", "missing");
        vars.setInfo("gtt1", gttype_string[gtt1]);
        is_match = false;
    }
    // only 2?
    else if (!is_call_1 && is_call_2)
    {
        vars.setInfo("type", "FP");
        vars.setInfo("kind", "missing");
        vars.setInfo("gtt2", gttype_string[gtt2]);
        is_match = false;
    }
    else if(is_call_1 && is_call_2)
    {
        // max number of alleles: 64 -- probably enough for diploid VCFs
        // (if there are more, preprocessing the VCFs/splitting out sample columns is
        //  probably the way to go)
        uint64_t alleles_seen_1 = 0;
        uint64_t alleles_seen_2 = 0;

        for (size_t gt = 0; gt < call1.ngt; ++gt)
        {
            uint64_t mask = 1;
            if (call1.gt[gt] > 0)
            {
                mask <<= (int8_t)call1.gt[gt];
            }
            alleles_seen_1 |= mask;
        }

        for (size_t gt = 0; gt < call2.ngt; ++gt)
        {
            uint64_t mask = 1;
            if (call2.gt[gt] > 0)
            {
                mask <<= (int8_t)call2.gt[gt];
            }
            alleles_seen_2 |= mask;
        }

        // see if we have the same alleles (over/undercall)
        std::string kind;
        if (alleles_seen_1 == alleles_seen_2)
        {
            if (gtt1 == gtt2)
            {
                // same alleles, same gt
                kind = "match";
                is_match = true;
            }
            else
            {
                kind = "gtmismatch";
                is_match = false;
            }
        }
        else
        {
            // mask out reference allele
            uint64_t nonref_als_1 = alleles_seen_1 & ((uint64_t)-2);
            uint64_t nonref_als_2 = alleles_seen_2 & ((uint64_t)-2);
            if (nonref_als_1 == nonref_als_2)
            {
                // different alleles, same gt
                kind = "gtmismatch";
                is_match = false;
            }
            else if ((nonref_als_1 & nonref_als_2) != 0)
            {
                kind = "alpartial";
                is_match = false;
            }
            else
            {
                kind = "almismatch";
                is_match = false;
            }
        }

        vars.setInfo("type", ( is_match ? "TP" : "FP" ));
        vars.setInfo("kind", kind.c_str());
        vars.setInfo("gtt1", gttype_string[gtt1]);
        vars.setInfo("gtt2", gttype_string[gtt2]);
    }
    else
    {
        // !is_call_1 && !is_call_2
        is_match = true;

        vars.setInfo("type", "N");
        vars.setInfo("kind", "match");
        vars.setInfo("gtt1", "no_call");
        vars.setInfo("gtt2", "no_call");
    }
#ifdef DEBUG_SIMPLECMP
    std::cerr << "                Calls match: " << is_match << "\n";
#endif

    if(is_match)
    {
        return dco_match;
    }
    else
    {
        return dco_mismatch;
    }
}

} // namespace haplotypes
