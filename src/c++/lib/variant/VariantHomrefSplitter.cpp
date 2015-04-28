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
 * \brief Homref interval splitter processor
 *
 * \file VariantHomrefSplitter.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantHomrefSplitter.hh"
#include "Error.hh"

#include <vector>
#include <list>

#define DEBUG_VARIANTHOMREFSPLITTER

namespace variant {


struct VariantHomrefSplitter::VariantHomrefSplitterImpl
{
    VariantHomrefSplitterImpl() {}

    VariantHomrefSplitterImpl(VariantHomrefSplitter::VariantHomrefSplitterImpl const & rhs)
    {
        buffered_variants = rhs.buffered_variants;
        output_variants = rhs.output_variants;
        vs = rhs.vs;
    }
    std::vector<Variants> buffered_variants;
    std::list<Variants> output_variants;
    Variants vs;

};


VariantHomrefSplitter::VariantHomrefSplitter()
{
    _impl = new VariantHomrefSplitterImpl();
}

VariantHomrefSplitter::VariantHomrefSplitter(VariantHomrefSplitter const & rhs)
{
    _impl = new VariantHomrefSplitterImpl(*rhs._impl);
}

VariantHomrefSplitter::~VariantHomrefSplitter()
{
    delete _impl;
}

VariantHomrefSplitter const & VariantHomrefSplitter::operator=(VariantHomrefSplitter const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    delete _impl;
    _impl = new VariantHomrefSplitterImpl(*rhs._impl);
    return *this;
}
        
/** enqueue a set of variants */
void VariantHomrefSplitter::add(Variants const & vs)
{
    if (vs.chr == _impl->vs.chr && vs.pos < _impl->vs.pos)
    {
        error("Variant added out of order at %s:%i / %i", vs.chr.c_str(), vs.pos, _impl->vs.pos);
    }

    _impl->buffered_variants.push_back(vs);
}

/**
 * @brief Return variant block at current position
 **/
Variants & VariantHomrefSplitter::current()
{
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantHomrefSplitter::advance()
{
    // refill output buffer?
    if (_impl->output_variants.empty()) {
        for (Variants & v : _impl->buffered_variants)
        {
#ifdef DEBUG_VARIANTHOMREFSPLITTER
            std::cerr << "VHRS: " << v << "\n";
#endif
            _impl->output_variants.push_back(v);
        }
        _impl->buffered_variants.clear();
    }

    if (_impl->output_variants.empty())
    {
        return false;
    }
    _impl->vs = _impl->output_variants.front();
    _impl->output_variants.pop_front();
    return true;
}

/** empty internal buffer */
void VariantHomrefSplitter::flush()
{
    _impl->buffered_variants.clear();
    _impl->output_variants.clear();
    _impl->vs = Variants();
}

}