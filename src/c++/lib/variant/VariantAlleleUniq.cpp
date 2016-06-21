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
 * Variant processing step to remove redundant alleles
 *
 *
 * \file VariantAlleleUniq.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantAlleleUniq.hh"

#include <map>

// #define DEBUG_VARIANTALLELEUNIQ

namespace variant
{

struct VariantAlleleUniq::VariantAlleleUniqImpl
{
    VariantAlleleUniqImpl() {}
    VariantAlleleUniqImpl(VariantAlleleUniqImpl const & rhs) : buffered_variants(rhs.buffered_variants), vs(rhs.vs) {}

    std::list<Variants> buffered_variants;
    Variants vs;
};

VariantAlleleUniq::VariantAlleleUniq()
{
    _impl = new VariantAlleleUniqImpl();
}

VariantAlleleUniq::VariantAlleleUniq(VariantAlleleUniq const & rhs) : _impl(new VariantAlleleUniqImpl(*rhs._impl))
{
}

VariantAlleleUniq::~VariantAlleleUniq()
{
    delete _impl;
}

VariantAlleleUniq const & VariantAlleleUniq::operator=(VariantAlleleUniq const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new VariantAlleleUniqImpl(*rhs._impl);
    return *this;
}

/** enqueue a set of variants */
void VariantAlleleUniq::add(Variants const & vs)
{
    if (vs.variation.size() <= 1 || vs.getInfoFlag("IMPORT_FAIL"))
    {
        _impl->buffered_variants.push_back(vs);
        return;
    }

    std::map<std::string, int > allele_map;
    std::vector<int> gt_mapping;
    gt_mapping.resize(vs.variation.size());
    int gt = 0;
    for (RefVar const & rv : vs.variation)
    {
        std::string rep = rv.repr();
        auto it = allele_map.find(rep);
        if (it == allele_map.end())
        {
            it = allele_map.emplace(rep, gt).first;
        }
        gt_mapping[gt] = it->second;
        ++gt;
    }

    if (allele_map.size() == vs.variation.size())
    {
        _impl->buffered_variants.push_back(vs);
        return;
    }

    Variants remapped = vs;
    remapped.variation.clear();
    gt = 0;
    std::vector<int> gt_mapping_2;
    int tmp = -1;
    gt_mapping_2.resize(vs.variation.size(), tmp);
    for (auto & p : allele_map)
    {
        remapped.variation.push_back(vs.variation[p.second]);
        gt_mapping_2[p.second] = gt++;
    }
    std::vector<int> gt_mapping_3;
    gt_mapping_3.resize(vs.variation.size());
    for (size_t i = 0; i < vs.variation.size(); ++i)
    {
        gt_mapping_3[i] = gt_mapping_2[gt_mapping[i]];
    }

#ifdef DEBUG_VARIANTALLELEUNIQ
    std::cerr << "GT remap at " << vs.chr << ":" << vs.pos << "\n";
    for (size_t i = 0; i < gt_mapping_3.size(); ++i)
    {
        std::cerr << i+1 << " -> " << (gt_mapping_3[i]+1) << "\n";
    }
#endif

    for (Call & c : remapped.calls)
    {
        for (size_t i = 0; i < c.ngt; ++i)
        {
            if (c.gt[i] > 0)
            {
                c.gt[i] = gt_mapping_3[c.gt[i]-1]+1;
            }
        }
    }

    for (auto & c : remapped.ambiguous_alleles)
    {
        for (auto & gt : c)
        {
            if (gt > 0)
            {
                gt = gt_mapping_3[gt-1];
            }
        }
    }
    _impl->buffered_variants.push_back(remapped);
}

/**
 * @brief Return variant block at current position
 **/
Variants & VariantAlleleUniq::current()
{
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantAlleleUniq::advance()
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
void VariantAlleleUniq::flush()
{
    _impl->buffered_variants.clear();
    _impl->vs = Variants();
}


}
