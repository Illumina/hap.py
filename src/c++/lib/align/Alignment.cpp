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
 *  \brief Alignment and Factory implementation
 *
 * \file Alignment.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Alignment.hh"

#include <cstdlib>
#include <cstring>
#include <sstream>

#include "Klib.hh"
#include "KlibGlobal.hh"
#include "helpers/Genetics.hh"

#include "Error.hh"

using namespace variant;
using namespace genetics;

Alignment::~Alignment() {}

/* get the best possible score for comparing two sequences of length len */
int Alignment::bestScore(int len)
{
    AlignmentParameters ap;
    getParameters(ap);
    return len*ap.maxScore();
}


Alignment * makeAlignment(const char * type)
{
    if(strstr(type, "klibg") == type)
    {
        return new KlibGlobalAlignment();
    }
    else if(strstr(type, "klib") == type)
    {
        return new KlibAlignment();
    }
    else
    {
        error("Unknown alignment type '%s'", type);
    }
    return NULL;
}

/**
 * @brief Format int encoded Cigar string
 *
 * @param tb for padding with "S" : begin
 * @param te for padding with "S" : end
 * @param altlen for padding with "S" : length of alternate sequence
 * @param ncigar length of cigar
 * @param cigar int* to cigar entries
 *
 */
std::string makeCigar(int tb, int te, int altlen, int ncigar, uint32_t * cigar)
{
    std::string cig;

    if (ncigar > 0)
    {
        std::ostringstream cigar_string;
        if (tb > 0)
        {
            cigar_string << tb << 'S';
        }

        for (int i = 0; i < ncigar; ++i)
        {
            cigar_string << (cigar[i] >> 4);
            uint8_t op = cigar[i] & 0x000f;
            switch(op)
            {
                case 0: cigar_string << 'M'; break;
                case 1: cigar_string << 'I'; break;
                case 2: cigar_string << 'D'; break;
            }
        }

        int end = altlen - te - 1;
        if (end > 0)
        {
            cigar_string << end << 'S';
        }

        cig = cigar_string.str();
    }

    return cig;
}


/** make variants from a cigar string */
void getVariantsFromCigar(std::string const & ref, std::string const & alt,
                          int r0, int a0,
                          uint32_t * cigar, int ncigar,
                          std::list<variant::RefVar> & target)
{

    bool have_rv = false;
    RefVar rv;
    rv.start = -1;
    rv.end = -1;
    rv.alt = "";

    int refpos = r0;
    int altpos = a0;

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = cigar[i] >> 4;
        uint8_t op = cigar[i] & 0x000f;
        switch(op)
        {
            case 0: // 'M'
            case 7: // '='
            case 8: // 'X'
                for(uint32_t j = 0; j < count; ++j)
                {
                    // check match / mismatch because ksw doesn't give this to us
                    if(ref[refpos] != alt[altpos])
                    {
                        // push a previous block which we can't append to
                        if(have_rv && rv.end < refpos-1)
                        {
                            target.push_back(rv);
                            have_rv = false;
                            rv.start = -1;
                            rv.end = -1;
                            rv.alt = "";
                        }

                        if(rv.start < 0)
                        {
                            rv.start = refpos;
                        }
                        rv.end = refpos;
                        rv.alt += alt[altpos];
                        have_rv = true;
                    }
                    else if(have_rv)
                    {
                        target.push_back(rv);
                        have_rv = false;
                        rv.start = -1;
                        rv.end = -1;
                        rv.alt = "";
                    }
                    ++refpos;
                    ++altpos;
                }
            break;
            case 1: // 'I' -> REF insertion in ALT = ALT deletion
                if(have_rv)
                {
                    target.push_back(rv);
                    have_rv = false;
                    rv.start = -1;
                    rv.end = -1;
                    rv.alt = "";
                }
                rv.start = refpos;
                rv.end = refpos + count - 1;
                rv.alt = "";
                have_rv = true;

                // shift the reference position
                refpos += count;
            break;
            case 2: // 'D' -> REF deletion = ALT insertion;
                if(have_rv)
                {
                    target.push_back(rv);
                    have_rv = false;
                    rv.start = -1;
                    rv.end = -1;
                    rv.alt = "";
                }
                // insert before reference pos

                // reference length = end - start + 1 == refpos-1 - refpos + 1 == 0
                // this is interpreted as an insertion before pos.
                rv.start = refpos;
                rv.end = refpos-1;
                rv.alt = alt.substr(altpos, count);

                have_rv = true;
                // shift the reference position up by one
                altpos += count;
            break;
        }
    }

    if(have_rv)
    {
        target.push_back(rv);
    }
}

/** get stats from a cigar string */
void getCigarStats(std::string const & ref, std::string const & alt,
                   int r0, int a0,
                   uint32_t * cigar, int ncigar,
                   int & softclipped, int & matches,
                   int & mismatches, int & ins, int & del)
{
    softclipped = r0;
    matches = 0;
    mismatches = 0;
    ins = 0;
    del = 0;
    int refpos = r0;
    int altpos = a0;

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = cigar[i] >> 4;
        uint8_t op = cigar[i] & 0x000f;
        switch(op)
        {
            case 0: // 'M'
            case 7: // '='
            case 8: // 'X'
                for(uint32_t j = 0; j < count; ++j)
                {
                    // check match / mismatch because ksw doesn't give this to us
                    if(ref[refpos] != alt[altpos])
                    {
                        ++mismatches;
                    }
                    else
                    {
                        ++matches;
                    }
                    ++refpos;
                    ++altpos;
                }
            break;
            case 1: // 'I' -> REF insertion in ALT = ALT deletion
                del += count;
                // shift the reference position
                refpos += count;
            break;
            case 2: // 'D' -> REF deletion = ALT insertion;
                ins += count;
                // shift the reference position up by one
                altpos += count;
            break;
        }
    }

    // TODO we might want to check if we're at the end of alt here.
    softclipped += ref.size() - refpos;
}

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param in_rv the RefVar record
 * @param aln the alignment interface to use
 * @param vars the primitive records
 */
void realignRefVar(FastaFile const & f, const char * chr, RefVar const & in_rv, Alignment * aln,
                   std::list<variant::RefVar> & vars)
{
    int64_t rstart = in_rv.start, rend = in_rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t)in_rv.alt.size();

    if(reflen < 2 || altlen < 2) //  || reflen == altlen)
    {
        // no complex ref / alt => use fast and simple function
        toPrimitives(f, chr, in_rv, vars);
        return;
    }

    std::string refseq = f.query(chr, rstart, rend);
    std::string altseq = in_rv.alt;

    aln->setRef(refseq.c_str());
    aln->setQuery(altseq.c_str());

    uint32_t * icigar;
    int r0, r1, a0, a1;
    int ncigar = 0;
    aln->getCigar(r0, r1, a0, a1, ncigar, icigar);

    RefVar rv;
    rv.start = -1;
    rv.end = -1;
    rv.alt = "";

    int refpos = r0;
    int altpos = a0;

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = icigar[i] >> 4;
        uint8_t op = (uint8_t) (icigar[i] & 0x000f);
        switch(op)
        {
            case 0: // 'M'
            case 7: // '='
            case 8: // 'X'
                for(uint32_t j = 0; j < count; ++j)
                {
                    // check match / mismatch because ksw doesn't give this to us
                    if(refseq[refpos] != altseq[altpos])
                    {
                        rv.start = rstart + refpos;
                        rv.end = rstart + refpos;
                        rv.alt = altseq[altpos];
                        vars.push_back(rv);
                    }
                    ++refpos;
                    ++altpos;
                }
            break;
            case 1: // 'I' -> REF insertion in ALT = ALT deletion
                rv.start = rstart + refpos;
                rv.end = rstart + refpos + count - 1;
                rv.alt = "";
                vars.push_back(rv);
                // shift the reference position
                refpos += count;
            break;
            case 2: // 'D' -> REF deletion = ALT insertion;
                // insert before reference pos

                // reference length = end - start + 1 == refpos-1 - refpos + 1 == 0
                // this is interpreted as an insertion before pos.
                rv.start = rstart + refpos;
                rv.end = rstart + refpos - 1;
                rv.alt = altseq.substr((unsigned long) altpos, count);
                vars.push_back(rv);

                altpos += count;
            break;
            default:break;
        }
    }
}

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param in_rv the RefVar record
 * @param snps the number of snps
 * @param ins the number of insertions
 * @param dels the number of deletions
 * @param homref the number of calls with no variation
 */
void realignRefVar(FastaFile const & f, const char * chr, variant::RefVar const & in_rv, Alignment * aln,
                   size_t & snps, size_t & ins, size_t & dels, size_t & homref,
                   size_t& transitions, size_t& transversions)
{
    int64_t rstart = in_rv.start, rend = in_rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t)in_rv.alt.size();

    if(reflen < 2 || altlen < 2)
    {
        // no complex ref / alt => use fast and simple function
        countRefVarPrimitives(f, chr, in_rv, snps, ins, dels, homref,
                              transitions, transversions);
        return;
    }

    std::string refseq = f.query(chr, rstart, rend);
    std::string altseq = in_rv.alt;

    aln->setRef(refseq.c_str());
    aln->setQuery(altseq.c_str());

    uint32_t * icigar;
    int r0, r1, a0, a1;
    int ncigar = 0;
    aln->getCigar(r0, r1, a0, a1, ncigar, icigar);

    int refpos = r0;
    int altpos = a0;
    bool isValidSnv(false);

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = icigar[i] >> 4;
        uint8_t op = (uint8_t) (icigar[i] & 0x000f);
        switch(op)
        {
            case 0: // 'M'
            case 7: // '='
            case 8: // 'X'
                for(uint32_t j = 0; j < count; ++j)
                {
                    // check match / mismatch because ksw doesn't give this to us
                    const char refBase(refseq[refpos]);
                    const char altBase(altseq[altpos]);

                    if (altBase != refBase)
                    {
                        ++snps;
                        const bool
                            isTransversion(snvIsTransversion(refBase, altBase,
                                                             isValidSnv));

                        if (isValidSnv) {
                            if (isTransversion) {
                                ++transversions;
                            } else {
                                ++transitions;
                            }
                        }
                    }
                    else
                    {
                        ++homref;
                    }
                    ++refpos;
                    ++altpos;
                }
            break;
            case 1: // 'I' -> REF insertion in ALT = ALT deletion
                dels += count;
                // shift the reference position
                refpos += count;
            break;
            case 2: // 'D' -> REF deletion = ALT insertion;
                // insert before reference pos

                // reference length = end - start + 1 == refpos-1 - refpos + 1 == 0
                // this is interpreted as an insertion before pos.
                ins+= count;
                altpos += count;
            break;
            default:break;
        }
    }
}

