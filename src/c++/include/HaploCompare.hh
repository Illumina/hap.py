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
 *  \brief Haplotype comparison definitions
 *
 *
 * \file HaploCompare.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Haplotype.hh"
#include "RefVar.hh"
#include "Alignment.hh"

#include <memory>
#include <vector>
#include <list>

namespace haplotypes
{

struct HaplotypeDiff
{
    HaplotypeDiff() : 
        score(-1), 
        hap1("."), hap2("."), 
        s1(-1), e1(-1), s2(-1), e2(-1), cigar("."),
        softclipped(0), matches(0), mismatches(0), ins(0), del(0) {}
    
    int score;

    std::string hap1;
    std::string hap2;

    int s1, e1, s2, e2;
    std::string cigar;

    // cigar stats
    int softclipped, matches, mismatches, ins, del;

    // variants
    std::list<variant::RefVar> vdiff;
};

std::ostream & operator<<(std::ostream & o, HaplotypeDiff const & oc);

/**
 * @brief Two-Haplotype Comparison
 */
struct HaploCompareImpl;
class HaploCompare
{
public:
    HaploCompare();
    ~HaploCompare();
    
    /**
     * Compare two haplotypes
     */

    /**
     * @brief Set to use a specific alignment object
     * 
     * The passed pointer will be deallocated when the HaploCompare object is destroyed.
     */
    void setAlignment(std::unique_ptr<Alignment> &&);
    Alignment* getAlignment();

    /**
     * Input sequences
     */
    void setRef(const char * );
    void setAlt(const char * );

    /**
     * Output 
     */

    /**
     * Use a global alignment to fully transform ref to alt.
     * 
     * If a local alignment will is used (a'la read realignment), 
     * returned variants will only be within the window r0 -> r1 and a0 -> a1
     * 
     * getPositions will return the locations of the fragments matched in ref and alt
     * (in case of a global alignment, these will be 0 and ref/alt length -1).
     * 
     * The score gives the alignment score for the local alignment within the positions.
     * 
     * Changes from ref to alt within the positions can be returned as a RefVar list.
     * 
     */
    void getVariants(std::list<variant::RefVar> & target);  

private:
    /** prevent assignment */
    virtual HaploCompare const & operator=(HaploCompare const &) { return *this; }

    HaploCompareImpl * _impl;
};

} // namespace haplotypes
