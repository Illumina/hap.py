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
 * \brief Haplotype comparison implementation
 *
 * \file HaploCompare.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "HaploCompare.hh"
#include "Alignment.hh"

#include <set>

using namespace variant;

namespace haplotypes
{

std::ostream & operator<<(std::ostream & o, HaplotypeDiff const & oc)
{
    o << oc.score << ":" 
      // << oc.hap1 << ":" 
      // << oc.hap2 << ":" 
      // << oc.cigar << ":" 
      << oc.softclipped << ":" 
      << oc.matches << ":" 
      << oc.mismatches << ":" 
      << oc.ins << ":" 
      << oc.del;

    size_t count = 0;
    for (RefVar const & rv : oc.vdiff)
    {
        if(count == 0)
        {
            o << ":";
        }
        o << rv.start << "-" << rv.end;
        if(rv.alt == "")
        {
            o << "=<DEL>";
        }
        else
        {
            o << "=" << rv.alt;
        }

        if(count + 2 <= oc.vdiff.size())
        {
            o << ",";
        }

        ++count;
    }
    return o;
}


struct HaploCompareImpl
{
    HaploCompareImpl() : align(makeAlignment("klibg")), valid_result(false)
    {}

    std::unique_ptr<Alignment> align;

    std::string ref;
    std::string alt;

    // output
    bool valid_result;
};

HaploCompare::HaploCompare()
{
    _impl = new HaploCompareImpl();
}

HaploCompare::~HaploCompare()
{
    delete _impl;
}
    
/**
 * @brief Set to use a specific alignment object
 * 
 * The passed pointer will be deallocated when the HaploCompare object is destroyed.
 */
void HaploCompare::setAlignment(std::unique_ptr<Alignment> && ap)
{
    _impl->align = std::move(ap);
}

Alignment* HaploCompare::getAlignment()
{
    return _impl->align.get();
}


/**
 * Input sequences
 */
void HaploCompare::setRef(const char * ref)
{
    _impl->align->setRef(ref);
    _impl->ref = ref;
}

void HaploCompare::setAlt(const char * query)
{
    _impl->align->setQuery(query);
    _impl->alt = query;
}

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
void HaploCompare::getVariants(std::list<RefVar> & target)
{
    if(_impl->ref.size() == 0 || _impl->alt.size() == 0)
    {
        return;
    }
    int r0, r1, a0, a1;
    uint32_t * cigar;
    int n_cigar;

    _impl->align->getCigar(r0, r1, a0, a1, n_cigar, cigar);

    getVariantsFromCigar(_impl->ref, _impl->alt, r0, a0, cigar, n_cigar, target);
}

} // namespace haplotypes
