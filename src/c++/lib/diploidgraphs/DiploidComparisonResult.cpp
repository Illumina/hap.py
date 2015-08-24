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
 * Comparison Result Implementation
 *
 * \file DiploidComparisonResult.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "DiploidComparisonResult.hh"

#include "HaploCompare.hh"
#include "DiploidReference.hh"
#include "Alignment.hh"

#include <memory>
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <sstream>
#include <iterator>

// #define DEBUG_DIPLOIDCOMPARISONRESULT

#ifdef DEBUG_DIPLOIDCOMPARISONRESULT
#include <fstream>
#endif

#include "Error.hh"

namespace haplotypes
{

std::ostream & operator<<(std::ostream & o, DiploidComparisonOutcome oc)
{
    switch(oc)
    {
        case dco_match:
            o << "dco_match";
            break;
        case dco_mismatch:
            o << "dco_mismatch";
            break;
        case dco_unknown:
        default:
            o << "dco_unknown";
            break;
    }
    return o;
}

std::ostream & operator<<(std::ostream & o, DiploidType oc)
{
    switch(oc)
    {
        case dt_hom:
            o << "hom";
            break;
        case dt_homref:
            o << "homref";
            break;
        case dt_het:
            o << "het";
            break;
        case dt_hetalt:
            o << "hetalt";
            break;
        default:
            o << "unknown";
            break;
    }
    return o;
}

std::ostream & operator<<(std::ostream & o, DiploidComparisonResult const & cr)
{
    o << cr.chr << "\t" << cr.start << "\t" << cr.end+1 << "\t" 
      << cr.outcome << "\t" 
      << cr.type1 << "\t" 
      << cr.type2 << "\t"
      << cr.n_paths1 << "\t" 
      << cr.n_paths2;
    
    if(cr.outcome == dco_mismatch)
    {
        o << "\t" << cr.n_pathsc;
        o << "\t" << cr.diffs[0];
        if(cr.type1 == dt_hetalt || cr.type2 == dt_hetalt)
        {
            o << "\t" << cr.diffs[1];
        }
        else
        {
            o << "\t.";
        }
    }
    else
    {
            o << "\t.\t.\t.";
    }

    return o;
}

void printDiploidComparisonResult(std::ostream & partial_bed, 
    DiploidComparisonResult const & dcr, bool output_sequences)
{
    partial_bed << 
        dcr.chr << "\t" << 
        dcr.start << "\t" <<
        dcr.end+1 << "\t" <<
        dcr.outcome << "\t" <<
        dcr.type1 << "\t" <<
        dcr.type2 << "\t" <<
        dcr.n_paths1 << ":" << dcr.n_paths2 << ":" << dcr.n_pathsc << ":" << dcr.n_nonsnp;
    
    if(dcr.outcome == dco_mismatch)
    {
        if(dcr.type1 == dt_hetalt || dcr.type2 == dt_hetalt)
        {
            int elen = dcr.end+1-dcr.start;
            if(elen <= 0)
            {
                elen = -1;
            }
            double identity_score = 0.5*dcr.diffs[0].score / double(elen) +
                                    0.5*dcr.diffs[1].score / double(elen);
            partial_bed << "\t" << identity_score
                        << "\t" << 2
                        << "\t" << (dcr.diffs[0].score + dcr.diffs[1].score);
            partial_bed << "\t" << dcr.diffs[0];
            if(output_sequences)
            {
                partial_bed << ":" << dcr.diffs[0].hap1
                            << ":" << dcr.diffs[0].hap2
                            << ":" << dcr.diffs[0].cigar;
            }
            partial_bed << "\t" << dcr.diffs[1];
            if(output_sequences)
            {
                partial_bed << ":" << dcr.diffs[1].hap1
                            << ":" << dcr.diffs[1].hap2
                            << ":" << dcr.diffs[1].cigar;
            }
        }
        else
        {
            int elen = dcr.end+1-dcr.start;
            if(elen <= 0)
            {
                elen = -1;
            }
            double identity_score = 0.5*dcr.diffs[0].score / double(elen);
            partial_bed << "\t" << identity_score
                        << "\t" << 1
                        << "\t" << dcr.diffs[0].score;
            partial_bed << "\t" << dcr.diffs[0];
            if(output_sequences)
            {
                partial_bed << ":" << dcr.diffs[0].hap1
                            << ":" << dcr.diffs[0].hap2
                            << ":" << dcr.diffs[0].cigar;
            }
        }
    }
    else if(dcr.outcome == dco_match)
    {
        partial_bed << "\t" << 1.0 <<   // combined identity score
                       "\t" << 1 <<     // number of haplotypes matched
                       "\t" << "." <<   // sum of alignment scores
                       "\t" << "." <<   // Haplotype 1
                       "\t" << ".";     // Haplotype 2
    }
    else
    {
        partial_bed << "\t" << "." <<   // combined identity score
                       "\t" << "." <<   // number of haplotypes matched
                       "\t" << "." <<   // sum of alignment scores
                       "\t" << "." <<   // Haplotype 1
                       "\t" << ".";     // Haplotype 2
    }
}

/**
 * Convert to JSON value
 */

Json::Value toJson(std::string const & s)
{
    return Json::Value(s);
}

Json::Value toJson(variant::RefVar const & rn)
{
    Json::Value val;
    val["pos"] = Json::Value::Int64(rn.start+1);
    val["end"] = Json::Value::Int64(rn.end+1);
    val["alt"] = rn.alt;
    val["flg"] = Json::Value::Int64(rn.flags);
    return val;
}


Json::Value toJson(AnnotatedRefVar const & rn)
{
    Json::Value val = toJson((variant::RefVar const &)rn);
    val["tags"] = toJson(rn.tags);
    return val;
}

Json::Value toJson(ReferenceNode const & rn)
{
    Json::Value val;
    val["pos"] = Json::Value::Int64(rn.start+1);
    val["end"] = Json::Value::Int64(rn.end+1);
    val["alt"] = rn.alt;
    val["clr"] = rn.color;
    val["flt"] = rn.filtered;
    val["het"] = rn.het;
    return val;
}

Json::Value toJson(HaplotypeDiff const & hd)
{
    Json::Value val;
    val["aln"] = hd.score;
    val["hs1"] = hd.hap1;
    val["hs2"] = hd.hap2;
    val["s1"] = hd.s1;
    val["s2"] = hd.s2;
    val["e1"] = hd.e1;
    val["e2"] = hd.e2;
    val["cig"] = hd.cigar;
    val["S"] = hd.softclipped;
    val["M"] = hd.matches;
    val["X"] = hd.mismatches;
    val["I"] = hd.ins;
    val["D"] = hd.del;

    if(hd.vdiff.size() > 0)
    {
        val["diff"] = toJson(hd.vdiff);
    }
    return val;
}

Json::Value toJson(DiploidComparisonResult const & cr)
{
    Json::Value val;

    val["chr"] = cr.chr;
    val["pos"] = Json::Value::Int64(cr.start+1);
    val["end"] = Json::Value::Int64(cr.end+1);
    val["ref"] = cr.refsq;

    {    
        std::ostringstream oss;
        oss << cr.outcome;
        val["res"] = oss.str();
    }
    {    
        std::ostringstream oss;
        oss << cr.type1;
        val["type1"] = oss.str();
    }
    {
        std::ostringstream oss;
        oss << cr.type2;
        val["type2"] = oss.str();
    }
    val["cmb"] = Json::Value();
    val["cmb"].append(Json::Value::Int64(cr.n_paths1));
    val["cmb"].append(Json::Value::Int64(cr.n_paths2));
    val["cmb"].append(Json::Value::Int64(cr.n_pathsc));
    val["nonsnp"] = Json::Value::Int64(cr.n_nonsnp);

    if(cr.type1 == dt_hetalt || cr.type2 == dt_hetalt)
    {
        val["hps"] = 2;
    }
    else
    {
        val["hps"] = 1;        
    }

    if(cr.outcome == dco_mismatch)
    {
        if(cr.type1 == dt_hetalt || cr.type2 == dt_hetalt)
        {
            val["aln"] = cr.diffs[0].score + cr.diffs[1].score;
            val["dfv"] = Json::Value();
            val["dfv"].append(toJson(cr.diffs[0]));
            val["dfv"].append(toJson(cr.diffs[1]));
        }
        else
        {
            val["aln"] = cr.diffs[0].score;
            val["dfv"] = Json::Value();
            val["dfv"].append(toJson(cr.diffs[0]));
        }
    }
    return val;
}

}
