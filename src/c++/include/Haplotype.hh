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
 *  \brief Haplotype enumeration interface
 *
 * \file Haplotype.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <string>
#include <list>

#include "RefVar.hh"

namespace haplotypes
{

/**
 * @brief Class to collect reference sequence modifications
 */
struct HaplotypeData;
class Haplotype
{
public:
    /**
     * @brief Initialize haplotype block for one chromosome
     * 
     * @param chr chromosome/contig to base the haplotype on
     * @param reference_file name of the reference sequence file
     */
    Haplotype(const char * chr, const char * reference_file);
    
    Haplotype(Haplotype const & );
    ~Haplotype();
    Haplotype const & operator=(Haplotype const &);

    /** add variant
     * 
     * Order of start coordinates must be increasing.
     */
    void addVar(int64_t start, int64_t end, std::string alt);
    void addVar(std::list<variant::RefVar> const & );

    /** Start and end of block */
    std::string chr() const;
    int64_t start() const;
    int64_t end() const;

    /**
     * Length difference compared to reference
     * 
     * min : minimum difference
     * max : maximum difference 
     *  (if multiple insertions/deletions are present)
     *  
     * sum : net difference in length
     * 
     */
    void lengthDiffs(int64_t & min, int64_t & max, int64_t & sum) const;

    /** Get the haplotype's modified reference sequence 
     * 
     * [-1, -1] : get only from start() to end()
     * 
     * Note that if start > start() or end < end(), they will be relative to 
     * the modified reference.
     * 
     */
    std::string seq(int64_t start = -1, int64_t end = -1) const;

    /**
     * @brief Canonical representation w.r.t. reference
     * 
     * Start and end can be used to give a region within which this block is canonical.
     */
    std::string repr(int64_t start = -1, int64_t end = -1) const;

    /**
     * @return true if no variants registered
     */
    bool noVar() const;

    /**
     * @brief Internal helper: clean up reference sequence hash (used in testing only)
     * 
     * To save time, we keep all reference_file Fastas in a hash. When feeding in temporary 
     * files, this is problematic since they can't be deleted until the Fasta object is 
     * deleted.
     */
    static void resetRefs();

private:
    HaplotypeData * _impl;
};

} // namespace haplotypes

