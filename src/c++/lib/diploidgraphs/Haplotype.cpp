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
 *  \brief Implementation of variant->Haplotype enumeration
 *
 * \file Haplotype.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Haplotype.hh"
#include "Error.hh"
#include "Fasta.hh"
#include "Alignment.hh"
#include "RefVar.hh"

#include <memory>
#include <list>
#include <map>
#include <sstream>
#include <cassert>
#include <limits>

using namespace variant;

namespace haplotypes
{


static std::map<std::string, std::shared_ptr<FastaFile> > FS_REF;

void Haplotype::resetRefs()
{
    FS_REF.clear();
}


/** Haplotype implementation data storage */
struct HaplotypeData
{
    HaplotypeData(std::string _chr, std::string refname) :
        chr(_chr), start(-1), end(-1)
    {
        auto rf = FS_REF.find(refname);
        if(rf == FS_REF.end())
        {
            std::shared_ptr<FastaFile> p(new FastaFile(refname.c_str()));
            rf = FS_REF.insert(std::make_pair(refname, p)).first;
        }
        refsq = rf->second;
    }

    HaplotypeData(HaplotypeData const & rhs)
    {
        *this = rhs;
    }
    HaplotypeData const & operator=(HaplotypeData const & rhs)
    {
        if(this == &rhs)
        {
            return *this;
        }
        refsq = rhs.refsq;
        chr = rhs.chr;
        start = rhs.start;
        end = rhs.end;
        v = rhs.v;
        return *this;
    }
    ~HaplotypeData() {}

    std::shared_ptr<FastaFile> refsq;
    std::string chr;
    int64_t start;
    int64_t end;

    std::list<RefVar> v;
};


Haplotype::Haplotype(const char * chr, const char * reference_file)
{
    _impl = new HaplotypeData(chr, reference_file);
}

Haplotype::Haplotype(Haplotype const & rhs) {
    _impl = new HaplotypeData(*rhs._impl);
}

Haplotype::~Haplotype()
{
    delete _impl;
}

Haplotype const & Haplotype::operator=(Haplotype const & rhs)
{
    *_impl = *(rhs._impl);
    return *this;
}

/** add variant */
void Haplotype::addVar(int64_t start, int64_t end, std::string alt)
{
    if(start <= _impl->end)
    {
        error("Haplotype variants out of order at %s:%i / %i", _impl->chr.c_str(), _impl->end, start);
    }
    int64_t reflen = end - start + 1;
    if(reflen < 0)
    {
        error("Reference length < 0 %s:%i / %i", _impl->chr.c_str(), start, end);
    }

    if(_impl->start < 0)
    {
        _impl->start = start;
    }

    _impl->end = end;
    
    RefVar rv;
    rv.start = start;
    rv.end = end;
    rv.alt = alt;
    _impl->v.push_back(rv);
}

void Haplotype::addVar(std::list<RefVar> const & l)
{
    for(RefVar const & rv : l)
    {
        addVar(rv.start, rv.end, rv.alt);
    }
}


/** Start and end of block */
std::string Haplotype::chr() const
{
    return _impl->chr;
}

int64_t Haplotype::start() const
{
    return _impl->start;
}

int64_t Haplotype::end() const
{
    return _impl->end;
}

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
void Haplotype::lengthDiffs(int64_t & min, int64_t & max, int64_t & sum) const
{
    min = std::numeric_limits<int64_t>::max();
    max = std::numeric_limits<int64_t>::min();
    sum = 0;
    int64_t shift = 0;
    for(RefVar const & rv : _impl->v)
    {
        int64_t reflen = rv.end - rv.start + 1;
        int64_t altlen = (int64_t)rv.alt.size();

        shift -= reflen - altlen;
        min = std::min(min, shift);
        max = std::max(max, shift);
    }    
    sum = shift;
}


/** get the full haplotype sequence */
std::string Haplotype::seq(int64_t start, int64_t end) const
{
    if (start < 0)
    {
        start = _impl->start;
    }
    if (end < 0)
    {
        end = _impl->end;
    }
    // is there no overlap?
    if(   start > std::max(_impl->start, _impl->end)
       || end < _impl->start || _impl->v.size() == 0)
    {
        // => return reference sequence
        return _impl->refsq->query(_impl->chr.c_str(), start, end);
    }

    // create modified reference

    std::string result;
    // reference length
    assert(_impl->end - _impl->start + 1 >= 0);
    int64_t istart = _impl->start;
    int64_t iend = _impl->end;
    if(_impl->end < _impl->start)
    {
        // this happens if we have only a single insertion
        result = _impl->refsq->query(_impl->chr.c_str(), istart, istart);
        iend = istart;
    }
    else
    {
        result = _impl->refsq->query(_impl->chr.c_str(), istart, iend);
    }
    int64_t shift = istart;
    for(RefVar const & rv : _impl->v)
    {
        int64_t pos = rv.start - shift;
        int64_t reflen = rv.end - rv.start + 1;
        int64_t altlen = (int64_t)rv.alt.size();

        assert(reflen >= 0);

        // support insertions with no reference character
        if(reflen == 0)
        {
            result.insert(pos, rv.alt);
        }
        else
        {
            result.replace(pos, reflen, rv.alt);
        }

        shift += reflen - altlen;
    }

    // these two cases are tricky because inside [istart -> iend], we don't actually have reference 
    // coordinates anymore. long deletions can remove basically all of the sequence at start
    if(start > istart)
    {   // overlapping, but starting after
        if (start >= istart + (signed)result.size())
        {
            result = "";
            istart = iend;
        }
        else
        {
            result = result.substr(start - istart);
            istart = start;
        }
    }
    if(end < iend)
    {
        if((signed)result.size() > end - istart + 1)
        {
            result = result.substr(0, (end - istart + 1));
        }
        iend = end;
    }

    if(end > iend)
    {
        result += _impl->refsq->query(_impl->chr.c_str(), iend+1, end);
    }
    
    if(start < istart)
    {   // overlapping but starting before
        result = _impl->refsq->query(_impl->chr.c_str(), start, istart-1) + result;
    }

    return result;
}

/**
 * @brief Canonical representation w.r.t. reference
 */
std::string Haplotype::repr(int64_t start, int64_t end) const
{
    std::ostringstream ss;
    if(_impl->start < 0 && _impl->end < 0)
    {
        // empty
        ss << _impl->chr << ":novar";
    }
    else
    {
        ss << _impl->chr << ":" << _impl->start << "-" << _impl->end << ":" << seq(start, end);
    }
    return ss.str();
}

/**
 * @return true if no variants registered
 */
bool Haplotype::noVar() const
{
    return (_impl->start < 0 && _impl->end < 0) || _impl->v.empty();
}


} // namespace haplotypes
