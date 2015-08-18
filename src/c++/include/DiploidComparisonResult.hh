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
 * Result data structures for diploid comparison
 *
 * \file DiploidComparisonResult.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#pragma once

#include <iostream>
#include "HaploCompare.hh"
#include "GraphReference.hh" 
#include "Variant.hh"

#include <json/json.h>

namespace haplotypes
{

typedef enum _DiploidComparisonOutcome
{
    dco_match,
    dco_mismatch,
    // failure when comparing
    dco_unknown
} DiploidComparisonOutcome;


std::ostream & operator<<(std::ostream & o, DiploidComparisonOutcome oc);

typedef enum _DiploidType
{
    dt_unknown,
    dt_het,
    dt_hetalt,
    dt_hom,
    dt_homref
} DiploidType;

/** construct DiploidType from het and homref fields in DiploidRef */
static inline DiploidType makeDiploidType(bool het, bool homref)
{
    if(het)
    {
        if(homref)
        {
            return dt_het;
        }
        else
        {
            return dt_hetalt;
        }
    }
    else
    {
        if(homref)
        {
            return dt_homref;
        }
        else
        {
            return dt_hom;
        }
    }    
}

std::ostream & operator<<(std::ostream & o, DiploidType oc);

struct AnnotatedRefVar : public variant::RefVar
{
    AnnotatedRefVar(variant::RefVar const & rv) : variant::RefVar(rv) {}

    std::list<std::string> tags;
};

struct DiploidComparisonResult
{
    // block position
    std::string chr;
    int64_t start;
    int64_t end;
    
    // outcome
    DiploidComparisonOutcome outcome;

    // combinatorics + performance: how many possibilities did we look at
    int64_t n_paths1, n_paths2, 
            n_pathsc; ///< and how many "expensive" comparisons / alignments did we do

    int64_t n_nonsnp;

    // matched haplotype pair types
    DiploidType type1, type2;

    // hap-cmp / DiploidCompare outputs
    // 
    // for matches, no diffs are given.
    // for hom<->hom(ref) mismatches, there will be one diff
    // for het<->hom(ref) mismatches, there will be one diff
    //    ... giving the diff of the alternate against the hom sequence
    // for het<->het mismatches, we give one diff 
    //    ... between the two alt calls
    //
    // If one of the inputs was het-alt, we need to give two diffs:
    //
    // for hetalt<->hom(ref) mismatches, there will be two diffs
    //    ... giving the diffs for each alternative against the hom sequence
    // for hetalt<->hetalt mismatches, there will be two diffs
    //    ... giving the diffs for the best two pairings.
    HaplotypeDiff diffs[2];

    std::string refsq;
};

/**
 * Bed-style output without sequences
 */
std::ostream & operator<<(std::ostream & o, DiploidComparisonResult const & cr);
void printDiploidComparisonResult(std::ostream & partial_bed, DiploidComparisonResult const & cr, bool output_sequences=true);

Json::Value toJson(std::string const & s);
Json::Value toJson(variant::RefVar const & rn);
Json::Value toJson(AnnotatedRefVar const & rn);
Json::Value toJson(ReferenceNode const & cr);
Json::Value toJson(HaplotypeDiff const & cr);
Json::Value toJson(DiploidComparisonResult const & cr);

/**
 * Convert to JSON values
 */
template <class _t1, class _t2>
static inline Json::Value toJson(std::pair<_t1, _t2> const & p)
{
	Json::Value v;
	v.append(toJson(p.first));
	v.append(toJson(p.second));
	return v;
}

template <class _t>
static inline Json::Value toJson(std::list<_t> const & l)
{
	Json::Value v;
	for(const _t & t: l)
	{
		v.append(toJson(t));
	}
	return v;
}

} // namespace haplotypes
