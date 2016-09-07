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
 * \brief Helper class to buffer RefVar objects into 'loci'
 *        (i.e. aggregate RefVar records by position)
 *
 * \file VariantAlleleSplitter.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#include <map>
#include <list>
#include <vector>

#include "variant/VariantAlleleSplitter.hh"
#include "Error.hh"

/* #define DEBUG_VARIANTALLELESPLITTER */

namespace variant {

struct VariantAlleleSplitter::VariantAlleleSplitterImpl
{
    VariantAlleleSplitterImpl() {}

    VariantAlleleSplitterImpl(VariantAlleleSplitter::VariantAlleleSplitterImpl const & rhs)
    {
        buffered_variants = rhs.buffered_variants;
        output_variants = rhs.output_variants;
        vs = rhs.vs;
    }
    std::vector<Variants> buffered_variants;
    std::list<Variants> output_variants;
    Variants vs;
};

VariantAlleleSplitter::VariantAlleleSplitter()
{
    _impl = new VariantAlleleSplitterImpl();
}

VariantAlleleSplitter::VariantAlleleSplitter(VariantAlleleSplitter const & rhs)
{
    _impl = new VariantAlleleSplitterImpl(*rhs._impl);
}

VariantAlleleSplitter::~VariantAlleleSplitter()
{
    delete _impl;
}

VariantAlleleSplitter const & VariantAlleleSplitter::operator=(VariantAlleleSplitter const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    delete _impl;
    _impl = new VariantAlleleSplitterImpl(*rhs._impl);
    return *this;
}

/** enqueue a set of variants */
void VariantAlleleSplitter::add(Variants const & vs)
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
 *
 * get_calls set to true to also retrieve GT information
 *
 * @param v Variant record to populate
 */
Variants & VariantAlleleSplitter::current()
{
#ifdef DEBUG_VARIANTALLELESPLITTER
    std::cerr << "Allele Splitter : Returning " << _impl->vs << "\n";
#endif
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantAlleleSplitter::advance()
{
    // refill output buffer?
    if (_impl->output_variants.empty()) {
        if (_impl->buffered_variants.empty())
        {
            return false;
        }
        size_t pos = 0;
        size_t n_samples = 0;
        std::string chr = "";

        struct CallInfo
        {
            CallInfo() : ad_ref(-1), ad_other(-1), dp(-1), nfilter(0) {}

            CallInfo(Call const & c) : ad_ref(c.ad_ref), ad_other(c.ad_other), dp(c.dp), qual(c.qual),
                                       nfilter(c.nfilter), formats(c.formats)
            {
                for (size_t i = 0; i < nfilter; ++i)
                {
                    filter[i] = c.filter[i];
                }
            }

            /** retain source records */
            int ad_ref, ad_other;
            int dp;
            float qual;
            size_t nfilter;
            std::string filter[MAX_FILTER];
            Json::Value formats;
        };

        // "half-call" -- one allele of a call
        struct HCall
        {
            RefVar rv;
            int ad;
            bool is_het;
            bool is_homref;
            size_t sample;
            Json::Value infos;
            CallInfo ci;
        };

        std::vector<HCall> loc_alleles;
        while(pos < _impl->buffered_variants.size())
        {
            Variants & v = _impl->buffered_variants[pos];
            chr = v.chr;

            n_samples = std::max(n_samples, v.calls.size());

            for (size_t i = 0; i < v.calls.size(); ++i)
            {
                Call & c = v.calls[i];
                if (c.ngt == 2 && c.gt[0] == c.gt[1]) // hom(alt)
                {
                    if (c.gt[0] > 0)
                    {
                        loc_alleles.push_back(HCall{v.variation[c.gt[0]-1], c.ad[0], false, false, i, v.infos, CallInfo(c)});
                    }
                    else if (c.gt[0] == 0)
                    {
                        RefVar hrv;
                        hrv.start = v.pos;
                        hrv.end = v.pos + v.len - 1;
                        hrv.alt = ".";
                        loc_alleles.push_back(HCall{hrv, c.ad_ref, false, true, i, v.infos, CallInfo(c)});
                    }
                }
                else // everything else
                {
                    bool any_nonref = false;
                    for (size_t j = 0; j < c.ngt; ++j)
                    {
                        if (c.gt[j] > 0)
                        {
                            loc_alleles.push_back(HCall{v.variation[c.gt[j]-1], c.ad[j], true, false, i, v.infos, CallInfo(c)});
                            any_nonref = true;
                        }
                    }

                    if(!any_nonref) // only add homref calls if there weren't any hets
                    {
                        for (size_t j = 0; j < c.ngt; ++j)
                        {
                            if (c.gt[j] == 0)
                            {
                                RefVar hrv;
                                hrv.start = v.pos;
                                hrv.end = v.pos + v.len - 1;
                                hrv.alt = ".";
                                loc_alleles.push_back(HCall{hrv, c.ad_ref, true, true, i, v.infos, CallInfo(c)});
                            }
                        }
                    }
                }
            }

            size_t sample = 0;
            for(auto li : v.ambiguous_alleles)
            {
                for(auto i : li)
                {
                    if(i > 0)
                    {
                        loc_alleles.push_back(HCall{v.variation[i-1], -1, true, false, sample, v.infos, CallInfo()});
                    }
                }
                ++sample;
            }

            ++pos;
        }
        if (pos > 0)
        {
            _impl->buffered_variants.erase(_impl->buffered_variants.begin(), _impl->buffered_variants.begin() + pos);
        }

        struct HCall_less
        {
            bool operator() (const HCall & l, const HCall & r)
            {
                bool result = false;
                if (l.rv.start < r.rv.start)
                {
                    result =  true;
                }
                else if (l.rv.start > r.rv.start)
                {
                    result =  false;
                }
                else if (l.rv.end > r.rv.end)
                {
                    result =  true;
                }
                else if (l.rv.end < r.rv.end)
                {
                    result =  false;
                }
                else
                {
                    result =  l.rv.alt < r.rv.alt;
                }
                return result;
            }
        };

        std::sort(loc_alleles.begin(), loc_alleles.end(), HCall_less());

#ifdef DEBUG_VARIANTALLELESPLITTER
        std::cerr << "Sorted loc_alleles in advance(): \n";
        for (auto const & l : loc_alleles)
        {
            std::cerr << "s" << l.sample << ": " << l.rv << "\n";
        }
#endif

        Variants cur;
        cur.chr = chr;
        cur.pos = -1;
        cur.len = 0;
        cur.variation.reserve(1);
        cur.calls.resize(n_samples);

        for (auto & p : loc_alleles)
        {
            // preserve INFO
            if (cur.pos < 0                                 // first one
             || (p.is_homref && cur.variation.size() > 0)   // var -> homref
             || (!p.is_homref && cur.variation.size() == 0) // new refvar
             || (cur.variation.size() > 0 && p.rv.repr() != cur.variation[0].repr())    // new refvar
             || cur.calls[p.sample].ngt > 0)                // call already there (!! multiple calls in same sample at same location)
            {
                if (cur.pos >= 0)
                {
                    _impl->output_variants.push_back(cur);
                }
                for (size_t i = 0; i < cur.calls.size(); ++i)
                {
                    cur.calls[i] = Call();
                }
                if (!p.is_homref)
                {
                    cur.variation.resize(1);
                    cur.variation[0] = p.rv;
                }
                else
                {
                    cur.variation.resize(0);
                }
                cur.pos = p.rv.start;
                cur.len = p.rv.end - p.rv.start + 1;
                cur.infos = p.infos;
            }
            for(auto const & m : p.infos.getMemberNames())
            {
                if(!cur.infos.isMember(m))
                {
                    cur.infos[m] = p.infos[m];
                }
            }
            cur.calls[p.sample].ngt = 2;
            cur.calls[p.sample].nfilter = p.ci.nfilter;
            cur.calls[p.sample].ad[0] = 0;
            cur.calls[p.sample].ad[1] = 0;
            cur.calls[p.sample].ad_ref = p.ci.ad_ref;
            cur.calls[p.sample].ad_other = p.ci.ad_other;
            cur.calls[p.sample].dp = p.ci.dp;
            cur.calls[p.sample].qual = p.ci.qual;
            cur.calls[p.sample].formats = p.ci.formats;
            for (size_t i = 0; i < p.ci.nfilter; ++i) cur.calls[p.sample].filter[i] = p.ci.filter[i];
            if(p.is_homref)
            {
                cur.calls[p.sample].gt[0] = 0;
                cur.calls[p.sample].gt[1] = 0;
                cur.calls[p.sample].ad[0] = p.ci.ad_ref;
                cur.calls[p.sample].ad[1] = p.ci.ad_ref;
            }
            else if (p.is_het)
            {
                cur.calls[p.sample].gt[0] = 0;
                cur.calls[p.sample].gt[1] = 1;
                cur.calls[p.sample].ad[0] = p.ci.ad_ref;
                cur.calls[p.sample].ad[1] = p.ad;
            }
            else
            {
                cur.calls[p.sample].gt[0] = 1;
                cur.calls[p.sample].gt[1] = 1;
                cur.calls[p.sample].ad[0] = p.ad;
                cur.calls[p.sample].ad[1] = p.ad;
            }
        }
        if (cur.pos >= 0)
        {
            _impl->output_variants.push_back(cur);
        }
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
void VariantAlleleSplitter::flush()
{
    _impl->buffered_variants.clear();
    _impl->output_variants.clear();
    _impl->vs = Variants();
}


}
