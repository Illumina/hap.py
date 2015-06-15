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
 * Partially phased diploid reference class. Enumerates possible haplotype pairs.
 *
 * \file DiploidReference.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Haplotype.hh"
#include "GraphReference.hh"
#include "Fasta.hh"

namespace haplotypes
{

/** Diploid (modified) reference sequence */
struct DiploidRef {
    bool het;
    bool homref;
    // if het==false, h2 is ignored
    std::string h1, h2;
    // reference sequence for this region
    std::string refsq;
};

std::ostream & operator<<(std::ostream & o, DiploidRef const & r);

struct DiploidReferenceImpl;
class DiploidReference
{
public:
    DiploidReference(GraphReference const & gr);
    DiploidReference(DiploidReference const &);
    DiploidReference const & operator=(DiploidReference const &);
    ~DiploidReference();

    /**
     * Get reference fasta file
     */
    FastaFile & getRefFasta();
    FastaFile const & getRefFasta() const;

    /**
     * Set the maximum number of paths to enumerate from the Graph reference
     */
    void setNPaths(int max_n_paths=-1);

    /**
     * Enumerate from set of Variants
     */
    void setRegion(
        const char * chr,
        int64_t start,
        int64_t end,
        std::list<variant::Variants> const & vars,
        int sample_ix = 0
    );

    /**
     * Enumeration: advance and return true if current is valid.
     *
     * Use like this:
     *
     * while(dr.hasNext())
     * {
     *     DiploidRef & hp (dr.next());
     *     // process hp
     *     ...
     *     dr.advance();
     * }
     *
     */
    bool hasNext();
    DiploidRef & next();
    void  advance();

    std::list<DiploidRef> const & result();

private:
	DiploidReferenceImpl * _impl;
};

}
