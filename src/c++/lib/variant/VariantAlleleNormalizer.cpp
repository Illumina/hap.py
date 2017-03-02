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
    VariantAlleleNormalizerImpl() : refpadding(false), homref(false), limit(-1), maxpos_chr(""), current_maxpos() {}
    VariantAlleleNormalizerImpl(VariantAlleleNormalizerImpl const & rhs) :
        buffered_variants(rhs.buffered_variants), reference(rhs.reference), refpadding(rhs.refpadding), homref(rhs.homref),
        limit(rhs.limit), maxpos_chr(rhs.maxpos_chr), current_maxpos(rhs.current_maxpos)
    {
        if (reference != "")
        {
            ref_fasta = std::move( std::unique_ptr<FastaFile>(new FastaFile( reference.c_str() ) ) );
        }
    }

    std::priority_queue<Variants, std::vector<Variants>, VariantCompare > buffered_variants;

    std::string reference;
    std::unique_ptr<FastaFile> ref_fasta;

    bool refpadding;
    bool homref;
    int64_t limit;
    std::string maxpos_chr;
    std::vector<int64_t> current_maxpos;

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
 * TODO : output homref as intervals
 */
int64_t VariantAlleleNormalizer::getLeftshiftLimit() const
{
    return _impl->limit;
}

void VariantAlleleNormalizer::setLeftshiftLimit(int64_t limit)
{
    _impl->limit = limit;
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

    // don't touch import fails
    if (vs.getInfoFlag("IMPORT_FAIL"))
    {
        _impl->buffered_variants.push(vs);
        return;
    }

    Variants nv(vs);
    size_t tmp = 0;
    _impl->current_maxpos.resize(std::max(_impl->current_maxpos.size(), nv.calls.size()), tmp);

#ifdef DEBUG_VARIANTNORMALIZER
    std::cerr << "before: " << nv << "\n";
#endif

    if(_impl->ref_fasta)
    {
        int64_t new_start = -1;
        int64_t new_end = -1;
        int64_t leftshift_limit = -1;

        if(_impl->limit >= 0)
        {
            leftshift_limit = nv.pos - _impl->limit;
        }

        for (size_t rvc = 0; rvc < nv.variation.size(); ++rvc)
        {
            int64_t this_leftshift_limit = leftshift_limit;
            if(nv.chr == _impl->maxpos_chr)
            {
                for(size_t j = 0; j < nv.calls.size(); ++j)
                {
                    for(size_t c = 0; c < nv.calls[j].ngt; ++c)
                    {
                        if(nv.calls[j].gt[c] - 1 == (int)rvc)
                        {
                            this_leftshift_limit = std::max(this_leftshift_limit, _impl->current_maxpos[j]);
                        }
                    }
                }
            }
#ifdef DEBUG_VARIANTNORMALIZER
            std::cerr << "leftshift limit for " << nv.variation[rvc] << " is " << this_leftshift_limit << "\n";
#endif
            leftShift(*(_impl->ref_fasta), nv.chr.c_str(), nv.variation[rvc], this_leftshift_limit, true);

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

        if (new_start >= 0 && new_end >= 0)
        {
            nv.pos = new_start;
            nv.len = new_end - new_start + 1;
            // handle insertions + corner case at start of chr
            if (nv.len == 0)
            {
                if(nv.pos > 0)
                {
                    --nv.pos;
                }
                nv.len = 1;
            }
        }
    }
#ifdef DEBUG_VARIANTNORMALIZER
    std::cerr << "after: " << nv << "\n";
#endif
    if(_impl->maxpos_chr != nv.chr)
    {
        for(size_t j = 0; j < nv.calls.size(); ++j)
        {
            if(!nv.calls[j].isNocall() && !nv.calls[j].isHomref())
            {
                _impl->current_maxpos[j] = nv.pos + nv.len - 1;
            }
        }
    }
    else
    {
        for(size_t j = 0; j < nv.calls.size(); ++j)
        {
            if(!nv.calls[j].isNocall() && !nv.calls[j].isHomref())
            {
                _impl->current_maxpos[j] = std::max(nv.pos + nv.len - 1, _impl->current_maxpos[j]);
            }
        }
    }
    _impl->maxpos_chr = nv.chr;
#ifdef DEBUG_VARIANTNORMALIZER
    std::cerr << "new max-shifting pos on " << _impl->maxpos_chr << " : ";
    for(size_t s = 0; s < _impl->current_maxpos.size(); ++s)
    {
        std::cerr << " s" << s << ": " << _impl->current_maxpos[s] << "  ";
    }
    std::cerr << "\n";
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
    _impl->buffered_variants = std::priority_queue<Variants, std::vector<Variants>, VariantCompare >();
    _impl->vs = Variants();
    _impl->maxpos_chr = "";
    _impl->current_maxpos.resize(0);
}

}
