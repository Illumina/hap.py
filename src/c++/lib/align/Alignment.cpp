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

#include "Error.hh"

using namespace variant;

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
 * @param n_cigar length of cigar
 * @param cigar int* to cigar entries
 * 
 */
std::string makeCigar(int tb, int te, int altlen, int n_cigar, uint32_t * cigar)
{
    std::string cig;

    if (n_cigar > 0) 
    {
        std::ostringstream cigar_string;
        if (tb > 0) 
        {
            cigar_string << tb << 'S';
        }

        for (int i = 0; i < n_cigar; ++i) 
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
                          uint32_t * cigar, int n_cigar,
                          std::list<variant::RefVar> & target)
{

    bool have_rv = false;
    RefVar rv;
    rv.start = -1;
    rv.end = -1;
    rv.alt = "";

    int refpos = r0;
    int altpos = a0;
    
    for (int i = 0; i < n_cigar; ++i)
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
                   uint32_t * cigar, int n_cigar,
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
    
    for (int i = 0; i < n_cigar; ++i)
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

