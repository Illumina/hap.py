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
 * \brief Diploid haplotype comparison implementation
 *
 * \file DiploidCompare.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "DiploidCompare.hh"
#include "HaploCompare.hh"
#include "DiploidReference.hh"
#include "Alignment.hh"

#include <memory>
#include <algorithm>
#include <map>
#include <set>
#include <limits>

#include "Error.hh"

//#define DEBUG_DIPLOIDCOMPARE

namespace haplotypes
{

struct DiploidCompareImpl
{
    DiploidCompareImpl(const char * ref_fasta) :
            dr(ref_fasta),
            nhap(4096),
            doAlignments(true)
    {
        matchScore = hcomp.getAlignment()->bestScore(1);
    }

    DiploidCompareImpl(DiploidCompareImpl const & rhs) :
            dr(rhs.dr),
            nhap(rhs.nhap),
            doAlignments(rhs.doAlignments)
    {
        matchScore = hcomp.getAlignment()->bestScore(1);
    }

    // diploid reference object
    DiploidReference dr;

    // parameters
    int nhap;

    HaploCompare hcomp;
    int matchScore;

    // updated every time we advance in nextResult
    DiploidComparisonResult cr;

    bool doAlignments;
};

DiploidCompare::DiploidCompare(const char * ref_fasta)
{
    _impl = new DiploidCompareImpl(ref_fasta);
}

DiploidCompare::DiploidCompare(DiploidCompare const & rhs)
{
    _impl = new DiploidCompareImpl(*rhs._impl);
}

DiploidCompare::~DiploidCompare()
{
    delete _impl;
}

DiploidCompare const & DiploidCompare::operator=(DiploidCompare const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new DiploidCompareImpl(*rhs._impl);
    return *this;
}

/**
 * Set the maximum number of haplotype blocks to enumerate before aborting
 */
void DiploidCompare::setMaxHapEnum(int nhap)
{
    _impl->nhap = nhap;
}

/**
 * Enable / disable the alignment step to find the best approximately matching
 * haplotypes (i.e. stop after match/mismatch status have been established)
 */
void DiploidCompare::setDoAlignments(bool doAlignments)
{
    _impl->doAlignments = doAlignments;
}

bool DiploidCompare::getDoAlignments()
{
    return _impl->doAlignments;
}

/**
 * @brief Set the region to compare in and reset the enumeration.
 *
 */
void DiploidCompare::setRegion(const char * chr, int64_t start, int64_t end,
                               std::list<variant::Variants> const & vars, int ix1, int ix2)
{
    _impl->cr.chr = chr;
    _impl->cr.start = start;
    _impl->cr.end = end;
    _impl->cr.refsq = ".";
    _impl->cr.outcome = dco_unknown;
    _impl->cr.type1 = dt_unknown;
    _impl->cr.type2 = dt_unknown;
    _impl->cr.n_paths1 = -1;
    _impl->cr.n_paths2 = -1;
    _impl->cr.n_pathsc = -1;
    _impl->cr.n_nonsnp = -1;
    _impl->cr.diffs[0] = HaplotypeDiff();
    _impl->cr.diffs[1] = HaplotypeDiff();

    // TODO: parameterize/fix
    if(end - start > 4096)
    {
        return;
    }
    _impl->dr.setNPaths(_impl->nhap);
    _impl->dr.setRegion(chr, start, end, vars, ix1);
    std::list<DiploidRef> di_haps1(_impl->dr.result());
    _impl->dr.setRegion(chr, start, end, vars, ix2);
    std::list<DiploidRef> di_haps2(_impl->dr.result());

#ifdef DEBUG_DIPLOIDCOMPARE
    std::cerr << "Input variants: " << "\n";
    for (auto const & v : vars) {
        std::cerr << v << "\n";
    }
    std::cerr << "At " << chr << ":" << start << "-" << end << "\n";
    std::cerr << "Pairs for truth: " << "\n";
    for(auto const & x : di_haps1)
    {
        std::cerr << x << "\n";
    }
    std::cerr << "Pairs for query: " << "\n";
    for(auto const & x : di_haps2)
    {
        std::cerr << x << "\n";
    }
#endif

    if(di_haps1.empty() || di_haps2.empty())
    {
        error("Unable to construct diploid haplotypes at %s:%i-%i", chr, start, end);
    }
    _impl->cr.n_paths1 = 2*di_haps1.size();
    _impl->cr.n_paths2 = 2*di_haps2.size();

    // this should be equal across all enumerated haps
    std::string refsq = di_haps1.front().refsq;
    _impl->cr.refsq = refsq;

    // this is where we store the haplotype sequences we have matched up
    // we will use this at the very end to find the variants they share
    // and the ones that are unique
    std::string matched_haplotypes_1[2];
    std::string matched_haplotypes_2[2];

    // 1. see if we can find a perfect allele match
    bool match_found = false;
    bool gt_mismatch_found = false; // detect cases where the alleles match, but the GT is wrong
    for (DiploidRef const & d1 : di_haps1)
    {
        DiploidType dt1 = makeDiploidType(d1.het, d1.homref);
        for (DiploidRef const & d2 : di_haps2)
        {
            DiploidType dt2 = makeDiploidType(d2.het, d2.homref);

            if(dt1 == dt2)
            {
                if(d1.het)  // d1.het == d2.het since dt1 == dt2
                {
                    bool direct_match = (d1.h1 == d2.h1 && d1.h2 == d2.h2);
                    bool flipped_match = (d1.h1 == d2.h2 && d1.h2 == d2.h1);
                    if( direct_match || flipped_match )
                    {
                        if(direct_match)
                        {
                            matched_haplotypes_1[0] = d1.h1;
                            matched_haplotypes_2[0] = d2.h1;
                            matched_haplotypes_1[1] = d1.h2;
                            matched_haplotypes_2[1] = d2.h2;
                        }
                        else
                        {
                            matched_haplotypes_1[0] = d1.h1;
                            matched_haplotypes_2[0] = d2.h2;
                            matched_haplotypes_1[1] = d1.h1;
                            matched_haplotypes_2[1] = d2.h2;
                        }
                        // het match
                        _impl->cr.type1 = dt1;
                        _impl->cr.type2 = dt2;
                        match_found = true;
                        break;
                    }
                }
                else
                {
                    // hom
                    if(d1.h1 == d2.h1)
                    {
                        matched_haplotypes_1[0] = d1.h1;
                        matched_haplotypes_2[0] = d2.h1;

                        _impl->cr.type1 = dt1;
                        _impl->cr.type2 = dt2;
                        match_found = true;
                        break;
                    }
                }
            }
            else if(dt1 != dt_hetalt && dt2 != dt_hetalt) // undercall / overcall?
            {
                std::string d1_alt;
                std::string d2_alt;

                if(d1.het)
                {
                    // one of the haplotypes in d1 is the reference sequence
                    if(d1.h1 == d1.refsq)
                    {
                        d1_alt = d1.h2;
                    }
                    else
                    {
                        d1_alt = d1.h1;
                    }
                }
                else
                {
                    d1_alt = d1.h1;
                }

                if(d2.het)
                {
                    // one of the haplotypes in d2 is the reference sequence
                    if(d2.h1 == d2.refsq)
                    {
                        d2_alt = d2.h2;
                    }
                    else
                    {
                        d2_alt = d2.h1;
                    }
                }
                else
                {
                    d2_alt = d2.h1;
                }

                if(d1_alt == d2_alt)
                {
                    matched_haplotypes_1[0] = d1_alt;
                    matched_haplotypes_2[0] = d2_alt;
                    gt_mismatch_found = true;
                    _impl->cr.type1 = dt1;
                    _impl->cr.type2 = dt2;
                    break;
                }
            }
            // hetalt vs het or hom doesn't need handling here
        }

        if(match_found)
        {
            _impl->cr.outcome = dco_match;
            _impl->cr.n_pathsc = 0;
            break;
        }
        if(gt_mismatch_found)
        {
            _impl->cr.outcome = dco_mismatch;
            _impl->cr.n_pathsc = 0;
            break;
        }
    }

    if(!_impl->doAlignments && !match_found && !gt_mismatch_found)
    {
        _impl->cr.outcome = dco_mismatch;
        return;
    }

    // no match was found: need to do in-depth comparison
    struct TmpResult
    {
        int score;
        int aln_count;
        DiploidType dt1, dt2;
        AlignmentResult alns[2];
        std::string haps1[2];
        std::string haps2[2];
    };

    if(!match_found && !gt_mismatch_found)
    {
        // find best match
        TmpResult best;
        best.score = std::numeric_limits<int>::min();
        best.aln_count = 0;

        for (DiploidRef const & d1 : di_haps1)
        {
            DiploidType dt1 = makeDiploidType(d1.het, d1.homref);
            for (DiploidRef const & d2 : di_haps2)
            {
                DiploidType dt2 = makeDiploidType(d2.het, d2.homref);

                if(dt1 != dt_hetalt && dt2 != dt_hetalt) // case 1: one alignment
                {
                    std::string d1_alt;
                    std::string d2_alt;

                    if(d1.het)
                    {
                        // one of the haplotypes in d1 is the reference sequence
                        if(d1.h1 == d1.refsq)
                        {
                            d1_alt = d1.h2;
                        }
                        else
                        {
                            d1_alt = d1.h1;
                        }
                    }
                    else
                    {
                        d1_alt = d1.h1;
                    }

                    if(d2.het)
                    {
                        // one of the haplotypes in d2 is the reference sequence
                        if(d2.h1 == d2.refsq)
                        {
                            d2_alt = d2.h2;
                        }
                        else
                        {
                            d2_alt = d2.h1;
                        }
                    }
                    else
                    {
                        d2_alt = d2.h1;
                    }

                    // we know d1.h1 != d2.h1 from above
                    // we do an alignment here, this is expensive, so we count how often this is done
                    ++_impl->cr.n_pathsc;
                    _impl->hcomp.setRef(d1_alt.c_str());
                    _impl->hcomp.setAlt(d2_alt.c_str());
                    int score = _impl->hcomp.getAlignment()->getScore();

                    // better than best known score?
                    // or less alignments
                    if(best.aln_count == 0 || best.aln_count > 1 || score > best.score)
                    {
                        best.aln_count = 1;
                        best.dt1 = dt1;
                        best.dt2 = dt2;
                        best.alns[0] = AlignmentResult(_impl->hcomp.getAlignment());
                        best.haps1[0] = d1_alt;
                        best.haps2[0] = d2_alt;
                    }
                }
                else if(best.aln_count != 1) // case 2: two alignments,
                                             // and we haven't found a single-pair alignment yet
                {
                    std::string d1_alt[2];
                    std::string d2_alt[2];

                    d1_alt[0] = d1.h1;
                    if(d1.het)
                    {
                        d1_alt[1] = d1.h2;
                    }
                    else
                    {
                        d1_alt[1] = d1.h1;
                    }

                    d2_alt[0] = d2.h1;
                    if(d2.het)
                    {
                        d2_alt[1] = d2.h2;
                    }
                    else
                    {
                        d2_alt[1] = d2.h1;
                    }

                    // compare pairs avoiding duplicate matches to
                    // one element
                    static const int ij[][4] = {
                        {0, 0, 1, 1},
                        {0, 1, 1, 0}
                    };

                    for (auto & ixs : ij)
                    {
                        int score = 0;
                        AlignmentResult ar[2];
                        for(int aln = 0; aln < 2; ++ aln)
                        {
                            int i = ixs[0 + 2*aln];
                            int j = ixs[1 + 2*aln];

                            ++_impl->cr.n_pathsc;
                            _impl->hcomp.setRef(d1_alt[i].c_str());
                            _impl->hcomp.setAlt(d2_alt[j].c_str());
                            score += _impl->hcomp.getAlignment()->getScore();
                            ar[aln] = AlignmentResult(_impl->hcomp.getAlignment());
                        }

                        if(best.aln_count == 0 || score > best.score)
                        {
                            best.score = score;
                            best.aln_count = 2;
                            best.dt1 = dt1;
                            best.dt2 = dt2;
                            best.alns[0] = ar[0];
                            best.alns[1] = ar[1];
                            best.haps1[0] = d1_alt[ixs[0]];
                            best.haps1[1] = d1_alt[ixs[2]];
                            best.haps2[0] = d2_alt[ixs[1]];
                            best.haps2[1] = d2_alt[ixs[3]];
                        }
                    }
                }
            }
        }

        _impl->cr.outcome = dco_mismatch;
        _impl->cr.type1 = best.dt1;
        _impl->cr.type2 = best.dt2;

        for (int i = 0; i < best.aln_count; ++i)
        {
            HaplotypeDiff & hd(_impl->cr.diffs[i]);

            matched_haplotypes_1[i] = best.haps1[i];
            matched_haplotypes_2[i] = best.haps2[i];

            hd.score = best.alns[i].score;
            hd.hap1 = best.haps1[i];
            hd.hap2 = best.haps2[i];

            hd.s1 = best.alns[i].s1;
            hd.e1 = best.alns[i].e1;
            hd.s2 = best.alns[i].s2;
            hd.e2 = best.alns[i].e2;

            hd.cigar = makeCigar(hd.s2, hd.e2, (int)hd.hap2.size(), best.alns[i].n_cigar, best.alns[i].cigar);

            getCigarStats(hd.hap1,
                          hd.hap2,
                          hd.s1, hd.s2,
                          best.alns[i].cigar, best.alns[i].n_cigar,
                          hd.softclipped,
                          hd.matches,
                          hd.mismatches,
                          hd.ins, hd.del);


            getVariantsFromCigar(hd.hap1, hd.hap2, hd.s1, hd.s2,
                                 best.alns[i].cigar, best.alns[i].n_cigar,
                                 hd.vdiff);
        }
    }
}

/**
 * @brief return (incremental) comparison outcome
 */
DiploidComparisonResult const & DiploidCompare::getResult()
{
    return _impl->cr;
}

} // namespace haplotypes
