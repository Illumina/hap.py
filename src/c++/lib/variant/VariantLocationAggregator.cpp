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
 *
 * Variant processing step to aggregate variants with the same start
 * location into single variant blocks
 *
 * \file VariantLocationAggregator.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantLocationAggregator.hh"

#include <list>

// #define DEBUG_VARIANTLOCATIONAGGREGATOR

namespace variant
{

struct VariantLocationAggregator::VariantLocationAggregatorImpl
{
    VariantLocationAggregatorImpl() : aggregationtype(VariantLocationAggregator::aggregate_nocall) {}
    VariantLocationAggregatorImpl(VariantLocationAggregatorImpl const & rhs) :
        buffered_variants(rhs.buffered_variants), aggregationtype(rhs.aggregationtype), vs(rhs.vs) {}
    std::list<Variants> buffered_variants;
    VariantLocationAggregator::AggregationType aggregationtype;
    Variants vs;
};

VariantLocationAggregator::VariantLocationAggregator()
{
    _impl = new VariantLocationAggregatorImpl();
}

VariantLocationAggregator::VariantLocationAggregator(VariantLocationAggregator const & rhs)
{
    _impl = new VariantLocationAggregatorImpl(*rhs._impl);
}

VariantLocationAggregator::~VariantLocationAggregator()
{
    delete _impl;
}

VariantLocationAggregator const & VariantLocationAggregator::operator=(VariantLocationAggregator const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new VariantLocationAggregatorImpl(*rhs._impl);
    return *this;
}

/** decide how to aggregate */
VariantLocationAggregator::AggregationType VariantLocationAggregator::getAggregationType() const
{
    return _impl->aggregationtype;
}

void VariantLocationAggregator::setAggregationType(VariantLocationAggregator::AggregationType at)
{
    _impl->aggregationtype = at;
}

/** enqueue a set of variants */
void VariantLocationAggregator::add(Variants const & vs)
{
#ifdef DEBUG_VARIANTLOCATIONAGGREGATOR
    std::cerr << vs << " / ";
    if( _impl->buffered_variants.empty())
    {
        std::cerr << "null";
    }
    else
    {
        std::cerr << _impl->buffered_variants.back();
    }
    std::cerr << "\n";
#endif
    if (   _impl->buffered_variants.empty()
        || vs.chr != _impl->buffered_variants.back().chr
        || vs.pos != _impl->buffered_variants.back().pos
        || vs.getInfoFlag("IMPORT_FAIL")
        || (   (vs.anyHomref() || _impl->buffered_variants.back().anyHomref())
            && (vs.pos + vs.len != _impl->buffered_variants.back().pos + _impl->buffered_variants.back().len) )
       )
    {
        _impl->buffered_variants.push_back(vs);
        return;
    }

    Variants & back = _impl->buffered_variants.back();
    std::list<size_t> combine;
    size_t nsamples = std::max(vs.calls.size(), back.calls.size());

    // check if we need to merge across types
    bool vs_has_snps = false, vs_has_indels = false;
    bool back_has_snps = false, back_has_indels = false;

    for (size_t i = 0; i < nsamples; ++i)
    {
        for(size_t ci = 0; ci < vs.calls[i].ngt; ++ci)
        {
            if(vs.calls[i].gt[ci] > 0)
            {
                RefVar const & rv = vs.variation[vs.calls[i].gt[ci]-1];
                int64_t reflen = rv.end - rv.start + 1;
                int64_t altlen = (int64_t)rv.alt.size();
                if(reflen == 1 && altlen == 1)
                {
                    vs_has_snps = true;
                }
                else
                {
                    vs_has_indels = true;
                }
            }
        }
        for(size_t ci = 0; ci < back.calls[i].ngt; ++ci)
        {
            if(back.calls[i].gt[ci] > 0)
            {
                RefVar const & rv = back.variation[back.calls[i].gt[ci]-1];
                int64_t reflen = rv.end - rv.start + 1;
                int64_t altlen = (int64_t)rv.alt.size();
                if(reflen == 1 && altlen == 1)
                {
                    back_has_snps = true;
                }
                else
                {
                    back_has_indels = true;
                }
            }
        }
    }

    bool compatible_vartype = true;
    if(( vs_has_snps && back_has_indels ) || ( vs_has_indels && back_has_snps ))
    {
        compatible_vartype = false;
    }
    if ( _impl->aggregationtype == aggregate_across_types )
    {
        compatible_vartype = true;
    }


    if(!compatible_vartype)
    {
        _impl->buffered_variants.push_back(vs);
        return;
    }

    for (size_t i = 0; i < nsamples; ++i)
    {
        bool can_combine = false;
        bool vs_has_call = (i < vs.calls.size()) && (!vs.calls[i].isNocall());
        bool back_has_call = (i < back.calls.size()) && (!back.calls[i].isNocall());

        bool compatible_filters = true;

        if(vs_has_call && back_has_call)
        {
            std::set<std::string> vs_filters;
            std::set<std::string> back_filters;

            for(size_t f = 0; f < vs.calls[i].nfilter; ++f)
            {
                vs_filters.insert(vs.calls[i].filter[f]);
            }
            for(size_t f = 0; f < back.calls[i].nfilter; ++f)
            {
                back_filters.insert(back.calls[i].filter[f]);
            }
            if(vs_filters != back_filters)
            {
                compatible_filters = false;
            }
        }

        switch(_impl->aggregationtype)
        {
            case aggregate_nocall:
                can_combine = !(back_has_call &&  vs_has_call);
                break;
            case aggregate_hetalt:
                can_combine = (!(back_has_call &&  vs_has_call))
                           || ( back_has_call && vs_has_call && back.calls[i].isHet() && vs.calls[i].isHet() );
                can_combine = can_combine && compatible_vartype && compatible_filters;
                break;
            case aggregate_ambigous:
                can_combine = compatible_filters;
                break;
            default:
            break;
        }
#ifdef DEBUG_VARIANTLOCATIONAGGREGATOR
        std::cerr << "\t" <<
                " vs_has_call: " << vs_has_call <<
                " can_combine: " << can_combine << "\n";
#endif

        if (can_combine)
        {
            combine.push_back(i);
        }
    }

#ifdef DEBUG_VARIANTLOCATIONAGGREGATOR
    std::cerr << "\t";
    for (size_t c : combine)
    {
        std::cerr << c << " ";
    }
    std::cerr << "\n";
#endif
    if (combine.empty())
    {
        _impl->buffered_variants.push_back(vs);
        return;
    }

    // combine info fields
    for(auto const & mn : vs.infos.getMemberNames())
    {
        if(!back.infos.isMember(mn))
        {
            back.infos[mn] = vs.infos[mn];
        }
    }
    // simple case: can merge all calls from vs into back
    Variants remaining_calls = vs;
    for (size_t c : combine)
    {
        bool het_call_vs = vs.calls[c].isHet();

        if (het_call_vs)
        {
            if (vs.calls[c].gt[0] > 0)
            {
                addAlleleToVariant(back, (int)c, vs.variation[vs.calls[c].gt[0]-1], true, vs.calls[c].ad[0]);
            }
            else if (vs.calls[c].gt[1] > 0)
            {
                addAlleleToVariant(back, (int)c, vs.variation[vs.calls[c].gt[1]-1], true, vs.calls[c].ad[1]);
            }
            else
            {
                error("Invalid het call at %s:%i", vs.chr.c_str(), vs.pos);
            }
        }
        else if(vs.calls[c].isHomalt())
        {
            addAlleleToVariant(back, (int)c, vs.variation[vs.calls[c].gt[0]-1], false, vs.calls[c].ad[0]);
        }
        else if(vs.calls[c].isHomref())
        {
            addRefAlleleToVariant(back, (int)c, vs.pos, vs.pos + vs.len - 1, false, std::max(vs.calls[c].ad_ref, vs.calls[c].ad[0]));
        }
        else
        {
            for (size_t i = 0; i < vs.calls[c].ngt; ++i)
            {
                int gt = vs.calls[c].gt[i];
                if(gt > 0)
                {
                    addAlleleToVariant(back, (int)c, vs.variation[gt-1], true, vs.calls[c].ad[i]);
                }
                else if(gt == 0)
                {
                    addRefAlleleToVariant(back, (int)c, vs.pos, vs.pos + vs.len - 1, false, std::max(vs.calls[c].ad_ref, vs.calls[c].ad[i]));
                }
            }
        }

        // merge ambiguous alleles, too
        if (c < vs.ambiguous_alleles.size())
        {
            int ado = vs.calls[c].ad_other;
            for (auto gt : vs.ambiguous_alleles[c])
            {
                if(gt > 0)
                {
                    addAlleleToVariant(back, (int)c, vs.variation[gt-1], true, ado);
                }
                else if(gt == 0)
                {
                   addRefAlleleToVariant(back, (int)c, vs.pos, vs.pos + vs.len - 1, false, ado);
                }
                // only add other depth only once
                ado = 0;
            }
            remaining_calls.ambiguous_alleles[c].clear();
        }

        // max depth / gq / qual
#ifdef DEBUG_VARIANTLOCATIONAGGREGATOR
        std::cerr << "\t" << "AD_REF: " << vs.calls[c].ad_ref << " / " << back.calls[c].ad_ref <<  "\n";
#endif
        back.calls[c].ad_ref = std::max(vs.calls[c].ad_ref, back.calls[c].ad_ref);
        back.calls[c].ad_other = std::max(vs.calls[c].ad_other, back.calls[c].ad_other);
        back.calls[c].dp = std::max(vs.calls[c].dp, back.calls[c].dp);
        back.calls[c].qual = std::max(vs.calls[c].qual, back.calls[c].qual);

        // have merged this one -> remove
        remaining_calls.calls[c].ngt = 0;
    }
#ifdef DEBUG_VARIANTLOCATIONAGGREGATOR
    std::cerr << "\t" << "Result: " << back << " / remaining: " << remaining_calls << "\n";
#endif

    // any calls remaining => add also
    for (size_t i = 0; i < remaining_calls.calls.size(); ++i)
    {
        if (remaining_calls.calls[i].ngt > 0
         || (remaining_calls.ambiguous_alleles.size() > i && remaining_calls.ambiguous_alleles[i].size() > 0))
        {
            _impl->buffered_variants.push_back(remaining_calls);
            break;
        }
    }
}

/**
 * @brief Return variant block at current position
 **/
Variants & VariantLocationAggregator::current()
{
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantLocationAggregator::advance()
{
    if (_impl->buffered_variants.empty())
    {
        return false;
    }
    else
    {
        _impl->vs = _impl->buffered_variants.front();
        _impl->buffered_variants.pop_front();
        return true;
    }
}

/** empty internal buffer */
void VariantLocationAggregator::flush()
{
    _impl->buffered_variants.clear();
    _impl->vs = Variants();
}



}
