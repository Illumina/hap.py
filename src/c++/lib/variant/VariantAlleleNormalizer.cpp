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
 * Variant processing step to normalize alleles by shifting and trimming
 *
 * \file VariantAlleleNormalizer.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantAlleleNormalizer.hh"
#include "Fasta.hh"

#include <queue>
#include <memory>

// #define DEBUG_VARIANTNORMALIZER

namespace variant
{

struct VariantAlleleNormalizer::VariantAlleleNormalizerImpl
{
    VariantAlleleNormalizerImpl() : refpadding(false), homref(false) {}
    VariantAlleleNormalizerImpl(VariantAlleleNormalizerImpl const & rhs) : 
        buffered_variants(rhs.buffered_variants), reference(rhs.reference), refpadding(rhs.refpadding), homref(rhs.homref)
    {
        if (reference != "")
        {
            ref_fasta = std::move( std::unique_ptr<FastaFile>(new FastaFile( reference.c_str() ) ) );
        }
    }

    struct Variants_less
    {
        bool operator() (Variants const & l, Variants const & r) 
        {
            if(l.pos != r.pos)
            {
                // std::priority_queue is a max-heap, we
                // want items with lowest pos to come out first
                return l.pos > r.pos;
            }
            else if (l.len != r.len)
            {
                // shorter ref allele comes out first
                return l.len > r.len;
            }
            else if(l.variation.size() != 1 || r.variation.size() != 1)
            {
                // if we don't have exactly one allele sort by number of alleles
                return l.variation.size() > r.variation.size();
            }
            else
            {
                // each one has exactly one allele, sort by this allele's alt
                return l.variation[0].alt > r.variation[0].alt;
            }
        }
    };

    std::priority_queue<Variants, std::vector<Variants>, Variants_less > buffered_variants;

    std::string reference;
    std::unique_ptr<FastaFile> ref_fasta;

    bool refpadding;
    bool homref;

    Variants vs;
};

VariantAlleleNormalizer::VariantAlleleNormalizer()
{
    _impl = new VariantAlleleNormalizerImpl();
}

VariantAlleleNormalizer::VariantAlleleNormalizer(VariantAlleleNormalizer const & rhs)
{
    _impl = new VariantAlleleNormalizerImpl(*rhs._impl);
}

VariantAlleleNormalizer::~VariantAlleleNormalizer()
{
    delete _impl;
}

VariantAlleleNormalizer const & VariantAlleleNormalizer::operator=(VariantAlleleNormalizer const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new VariantAlleleNormalizerImpl(*rhs._impl);
    return *this;
}
 
void VariantAlleleNormalizer::setReference(std::string const & fasta)
{
    _impl->reference = fasta;
    _impl->ref_fasta = std::move( std::unique_ptr<FastaFile>(new FastaFile(_impl->reference.c_str() ) ) );  
}


/**
 * @brief Enable / disable reference padding
 * 
 * RefVars can be trimmed to include / not include reference bases
 */
bool VariantAlleleNormalizer::getEnableRefPadding() const
{
    return _impl->refpadding;
}

void VariantAlleleNormalizer::setEnableRefPadding(bool padding)
{
    _impl->refpadding = padding;
}

/**
 * @brief Enable / disable returning of homref intervals
 * 
 */
bool VariantAlleleNormalizer::getEnableHomrefVariants() const
{
    return _impl->homref;
}

void VariantAlleleNormalizer::setEnableHomrefVariants(bool homref)
{
    _impl->homref = homref;
}

/** enqueue a set of variants */
void VariantAlleleNormalizer::add(Variants const & vs)
{
    bool all_homref = true;
    for(Call const & c : vs.calls)
    {
        if(!(c.isHomref() || c.isNocall()))
        {
            all_homref = false;
            break;
        }
    }

    for(auto const & x : vs.ambiguous_alleles)
    {
        if(!x.empty())
        {
            all_homref = false;
            break;                    
        }
    }

    if (all_homref && !_impl->homref)
    {
#ifdef DEBUG_VARIANTNORMALIZER
        std::cerr << "Skipping homref variant: " << vs << "\n";
#endif
        return;
    }

    Variants nv(vs);

#ifdef DEBUG_VARIANTNORMALIZER
    std::cerr << "before: " << nv << "\n";
#endif
    if (_impl->ref_fasta)
    {
        int64_t new_start = -1;
        int64_t new_end = -1;

        for (size_t rvc = 0; rvc < nv.variation.size(); ++rvc)
        {
            leftShift(*(_impl->ref_fasta), nv.chr.c_str(), nv.variation[rvc]);
            trimLeft(*(_impl->ref_fasta), nv.chr.c_str(), nv.variation[rvc], _impl->refpadding);
            trimRight(*(_impl->ref_fasta), nv.chr.c_str(), nv.variation[rvc], _impl->refpadding);            
            if(new_start < 0 || new_start > nv.variation[rvc].start)
            {
                new_start = nv.variation[rvc].start;
            }
            if(new_end < 0 || new_end > nv.variation[rvc].end)
            {
                new_end = nv.variation[rvc].end;
            }
        }

        if (new_start > 0 && new_end > 0)
        {
            nv.pos = new_start;
            nv.len = new_end - new_start + 1;
        }
    }
#ifdef DEBUG_VARIANTNORMALIZER
    std::cerr << "after: " << nv << "\n";
#endif
    _impl->buffered_variants.push(nv);
}

/**
 * @brief Return variant block at current position
 **/
Variants & VariantAlleleNormalizer::current()
{
#ifdef DEBUG_VARIANTNORMALIZER
    std::cerr << "returning " << _impl->vs << "\n";
#endif
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantAlleleNormalizer::advance()
{
    if (_impl->buffered_variants.empty())
    {
#ifdef DEBUG_VARIANTNORMALIZER
        std::cerr << "no buffered variants\n";
#endif
        return false;
    }
    else
    {
        _impl->vs = _impl->buffered_variants.top();
        _impl->buffered_variants.pop();
#ifdef DEBUG_VARIANTNORMALIZER
        std::cerr << "Variants left: " << _impl->buffered_variants.size() << " / empty: " << _impl->buffered_variants.empty() <<  "\n";
#endif
        return true;        
    }
}

/** empty internal buffer */
void VariantAlleleNormalizer::flush()
{
    _impl->buffered_variants = std::priority_queue<Variants, std::vector<Variants>, VariantAlleleNormalizerImpl::Variants_less >();
    _impl->vs = Variants();
}

}
