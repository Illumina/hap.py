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

#include <queue>
#include <vector>
#include <list>

// #define DEBUG_VARIANTHOMREFSPLITTER

namespace variant {

namespace ns_varianthomrefsplitter
{
    typedef std::priority_queue<
        Variants,
        std::vector<Variants>,
        VariantCompare
    > VariantQueue;
}

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

    ns_varianthomrefsplitter::VariantQueue output_variants;

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
    if (_impl->buffered_variants.size() > 0 &&
        vs.chr == _impl->buffered_variants.back().chr &&
        vs.pos < _impl->buffered_variants.back().pos)
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
#ifdef DEBUG_VARIANTHOMREFSPLITTER
    std::cerr << "VHRS-output: " << _impl->vs << "\n";
#endif
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantHomrefSplitter::advance()
{
    // refill output buffer?
    if (_impl->output_variants.empty())
    {
        for (Variants & v : _impl->buffered_variants)
        {
            if (v.anyHomref())
            {
                Variants non_hr = v;
                int n_non_hr = v.calls.size();
                for (size_t q = 0; q < v.calls.size(); ++q)
                {
                    if(v.calls[q].isHomref())
                    {
                        non_hr.calls[q] = Call();
                        --n_non_hr;
                    }
                    else if(v.calls[q].isNocall())
                    {
                        --n_non_hr;
                    }
                    else
                    {
                        v.calls[q] = Call();
                    }
                }
                if (n_non_hr || non_hr.anyAmbiguous())
                {
#ifdef DEBUG_VARIANTHOMREFSPLITTER
                    std::cerr << "VHRS-non-hr-pass-on: " << v << "\n";
#endif
                    _impl->output_variants.push(non_hr);
                }

                v.variation.clear();
                for(auto & l : v.ambiguous_alleles)
                {
                    l.clear();
                }

#ifdef DEBUG_VARIANTHOMREFSPLITTER
                std::cerr << "VHRS-hr-split: " << v << "\n";
#endif
                // TODO handle very long blocks?
                int64_t end_pos = v.pos + v.len;
                v.len = 1;
                while(v.pos < end_pos)
                {
                    _impl->output_variants.push(v);
                    ++v.pos;
                }
            }
            else
            {
#ifdef DEBUG_VARIANTHOMREFSPLITTER
                std::cerr << "VHRS-pass-on: " << v << "\n";
#endif
                _impl->output_variants.push(v);
            }
        }
        _impl->buffered_variants.clear();
    }

    if (_impl->output_variants.empty())
    {
        return false;
    }
    _impl->vs = _impl->output_variants.top();
    _impl->output_variants.pop();
    return true;
}

/** empty internal buffer */
void VariantHomrefSplitter::flush()
{
    _impl->buffered_variants.clear();
    _impl->output_variants = ns_varianthomrefsplitter::VariantQueue();
    _impl->vs = Variants();
}

}
