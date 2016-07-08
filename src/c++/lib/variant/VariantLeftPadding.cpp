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
 * \file VariantLeftpadding.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantLeftpadding.hh"
#include "Error.hh"

#include <queue>
#include <vector>
#include <list>

// #define DEBUG_VARIANTLEFTPADDING

namespace variant {

namespace ns_variantleftpadding
{
    typedef std::priority_queue<
        Variants,
        std::vector<Variants>,
        VariantCompare
    > VariantQueue;
}

struct VariantLeftpadding::VariantLeftpaddingImpl
{
    VariantLeftpaddingImpl() {}

    VariantLeftpaddingImpl(VariantLeftpadding::VariantLeftpaddingImpl const & rhs)
    {
        buffered_variants = rhs.buffered_variants;
        output_variants = rhs.output_variants;
        ref = rhs.ref;
        vs = rhs.vs;
    }

    std::vector<Variants> buffered_variants;

    std::shared_ptr<FastaFile> ref;

    ns_variantleftpadding::VariantQueue output_variants;

    Variants vs;
};


VariantLeftpadding::VariantLeftpadding()
{
    _impl = new VariantLeftpaddingImpl();
}

VariantLeftpadding::VariantLeftpadding(VariantLeftpadding const & rhs)
{
    _impl = new VariantLeftpaddingImpl(*rhs._impl);
}

VariantLeftpadding::~VariantLeftpadding()
{
    delete _impl;
}

VariantLeftpadding const & VariantLeftpadding::operator=(VariantLeftpadding const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    delete _impl;
    _impl = new VariantLeftpaddingImpl(*rhs._impl);
    return *this;
}

/**
 * @brief Reference for shifting
 *
 */
void VariantLeftpadding::setReference(std::string const & fasta)
{
    ref = std::make_shared(new FastaFile(fasta.c_str()));
}

/** enqueue a set of variants */
void VariantLeftpadding::add(Variants const & vs)
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
Variants & VariantLeftpadding::current()
{
#ifdef DEBUG_VARIANTLEFTPADDING
    std::cerr << "VLP-output: " << _impl->vs << "\n";
#endif
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantLeftpadding::advance()
{
    // refill output buffer?
    if (_impl->output_variants.empty())
    {
        for (Variants & v : _impl->buffered_variants)
        {
            bool do_pad = false;
            for(auto & alt : v.variation)
            {
                const int64_t reflen = alt.end - alt.start;
                const int64_t altlen = alt.alt.size();

                if(reflen == 0 || altlen == 0)
                {
                    do_pad = true;
                    break;
                }
            }

            if(do_pad)
            {
                v.pos--;
                std::string padding = _impl->ref.query(v.chr, v.pos, v.pos);
                for(auto & al : v.variation)
                {
                    al.start--;
                }
            }

            _impl->output_variants.push(v);
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
void VariantLeftpadding::flush()
{
    _impl->buffered_variants.clear();
    _impl->output_variants = ns_varianthomrefsplitter::VariantQueue();
    _impl->vs = Variants();
}

}
