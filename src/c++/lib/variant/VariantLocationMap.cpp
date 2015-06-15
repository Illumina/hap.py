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
 * Functions for aggregating variants by location
 *
 * \file VariantLocationMap.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Variant.hh"


namespace variant
{

/**
 * @brief Combine alleles in a single Variants record
 *
 */
void addAlleleToVariant(Variants & target, int sample, RefVar const & rv, bool is_het, int ad)
{
    int64_t new_start = std::min(rv.start, target.pos);
    int64_t new_end = std::max(rv.end, target.pos + target.len - 1);
    target.pos = new_start;
    target.len = new_end - new_start + 1;
    target.variation.push_back(rv);
    size_t gt_id = target.variation.size();

    if ((int)target.calls.size() <= sample)
    {
        target.calls.resize(sample+1);
    }

    if ((int)target.ambiguous_alleles.size() <= sample)
    {
        target.ambiguous_alleles.resize(sample+1);
    }

    if (target.ambiguous_alleles[sample].size() > 0)
    {
        target.ambiguous_alleles[sample].push_back(gt_id);
        if (!is_het)
        {
            target.ambiguous_alleles[sample].push_back(gt_id);
        }
        target.calls[sample].ad_other += ad;
    }
    else
    {
        Call & call = target.calls[sample];

        if (is_het)
        {
            // het -> find last gtt that is zero
            int zix = -1;
            // at least one
            call.ngt = std::min(static_cast<size_t>(MAX_GT), std::max(static_cast<size_t>(2), call.ngt + 1));
            for (int i = 0; i < MAX_GT; ++i)
            {
                if (call.gt[i] < 0)
                {
                    call.gt[i] = 0;
                }
                if(call.gt[i] == 0)
                {
                    zix = i;
                }
            }
            if (zix < 0)
            {
#ifdef VARIANTLOCATIONMAP_WARN_AMBIGUOUS
                fprintf(stderr, "[W] Too many het alleles at %s:%i (%i)\n",
                        target.chr.c_str(), (int)target.pos, (int)call.ngt + 1);
#endif
                call.ad_other = std::max(0, ad);
                for (size_t i = 0; i < call.ngt; ++i)
                {
                    if (call.gt[i] >= 0)
                    {
                        target.ambiguous_alleles[sample].push_back(call.gt[i]);
                        call.gt[i] = -1;
                        call.ad_other += std::max(0, call.ad[i]);
                        call.ad[i] = 0;
                    }
                }
                target.ambiguous_alleles[sample].push_back(gt_id);
            }
            else
            {
                call.gt[zix] = gt_id;
                call.ad[zix] = ad;
            }
        }
        else
        {
            // hom -> this needs to be the only refvar here
            if (call.ngt > 0)
            {
#ifdef VARIANTLOCATIONMAP_WARN_AMBIGUOUS
                fprintf(stderr, "[W] Too many hom alleles at %s:%i / %i (%i)\n",
                        target.chr.c_str(), (int)target.pos, (int)sample, (int)call.ngt + 1);
#endif
                call.ad_other = std::max(0, ad);
                for (size_t i = 0; i < call.ngt; ++i)
                {
                    if (call.gt[i] >= 0)
                    {
                        target.ambiguous_alleles[sample].push_back(call.gt[i]);
                        call.gt[i] = -1;
                        call.ad_other += std::max(0, call.ad[i]);
                        call.ad[i] = 0;
                    }
                }
                target.ambiguous_alleles[sample].push_back(gt_id);
                target.ambiguous_alleles[sample].push_back(gt_id);
            }
            else
            {
                call.ngt = 2;
                call.gt[0] = gt_id;
                call.gt[1] = gt_id;
                call.ad[0] = ad;
                call.ad[1] = ad;
            }
        }
    }
}

/**
 * @brief Add ref block to variant record
 */
void addRefAlleleToVariant(Variants & target, int sample, int64_t start, int64_t end, bool is_het, int ad)
{
    int64_t new_start = std::min(start, target.pos);
    int64_t new_end = std::max(end, target.pos + target.len - 1);
    target.pos = new_start;
    target.len = new_end - new_start + 1;

    if ((int)target.calls.size() <= sample)
    {
        target.calls.resize(sample+1);
    }
    if ((int)target.ambiguous_alleles.size() <= sample)
    {
        target.ambiguous_alleles.resize(sample+1);
    }

    if (target.ambiguous_alleles[sample].size() > 0)
    {
        target.ambiguous_alleles[sample].push_back(0);
        if (!is_het)
        {
            target.ambiguous_alleles[sample].push_back(0);
        }

        target.calls[sample].ad_other += ad;
    }
    else
    {
        Call & call = target.calls[sample];

        // max ref coverage
        call.ad_ref = std::max(call.ad_ref, ad);

        if (is_het)
        {
            // het -> find last gtt that is zero
            int zix = -1;
            // don't replace allele when call is already homref
            if (!(call.isHomref() && call.ngt == 2))
            {
                call.ngt = std::min(static_cast<size_t>(MAX_GT), std::max(static_cast<size_t>(2), call.ngt + 1));
                for (int i = 0; i < MAX_GT; ++i)
                {
                    if(call.gt[i] == 0)
                    {
                        zix = i;
                    }
                }
            }
            // if the call is already homref, don't add more ref alleles
            if (zix < 0)
            {
#ifdef VARIANTLOCATIONMAP_WARN_AMBIGUOUS
                fprintf(stderr, "[W] Too many het alleles at %s:%i (%i)\n",
                        target.chr.c_str(), (int)target.pos, (int)call.ngt + 1);
#endif
                call.ad_other = std::max(ad, 0);
                for (size_t i = 0; i < call.ngt; ++i)
                {
                    if (call.gt[i] >= 0)
                    {
                        target.ambiguous_alleles[sample].push_back(call.gt[i]);
                        call.gt[i] = -1;
                        call.ad_other += std::max(0, call.ad[i]);
                        call.ad[i] = 0;
                    }
                }
                target.ambiguous_alleles[sample].push_back(0);
            }
            else
            {
                call.gt[zix] = 0;
                call.ad[zix] = ad;
            }
        }
        else
        {
            // hom -> this needs to be the only refvar here
            if (call.ngt > 0)
            {
#ifdef VARIANTLOCATIONMAP_WARN_AMBIGUOUS
                fprintf(stderr, "[W] Too many hom alleles at %s:%i / %i (%i)\n",
                        target.chr.c_str(), (int)target.pos, (int)sample, (int)call.ngt + 1);
#endif
                call.ad_other = std::max(ad, 0);
                for (size_t i = 0; i < call.ngt; ++i)
                {
                    if (call.gt[i] >= 0)
                    {
                        target.ambiguous_alleles[sample].push_back(call.gt[i]);
                        call.gt[i] = -1;
                        call.ad_other += std::max(0, call.ad[i]);
                        call.ad[i] = 0;
                    }
                }
                target.ambiguous_alleles[sample].push_back(0);
                target.ambiguous_alleles[sample].push_back(0);
            }
            else
            {
                // at least one
                call.ngt = 2;
                call.gt[0] = 0;
                call.gt[1] = 0;
                call.ad[0] = ad;
                call.ad[1] = ad;
            }
        }
    }
}

VariantLocationMap::iterator addToLocationMap(
    VariantLocationMap & locmap, int sample,
    RefVar const & rv, bool is_het)
{
    auto cit = locmap.find(rv.start);
    if (cit == locmap.end())
    {
        cit = locmap.emplace(rv.start, Variants()).first;
        cit->second.pos = rv.start;
        cit->second.len = rv.end - rv.start + 1;
    }
    else
    {
        // cit->pos should be equal to rv.start because we index by this in the map
        // cit->pos = rv.start;
        cit->second.len = std::max(rv.end - rv.start + 1, cit->second.len);
    }

    addAlleleToVariant(cit->second, sample, rv, is_het);

    return cit;
}

}
