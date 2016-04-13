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
 * Track named regions
 *
 * \file QuantifyRegions.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_QUANTIFYREGIONS_HH
#define HAPLOTYPES_QUANTIFYREGIONS_HH

#include <memory>
#include <map>
#include <string>

#include <htslib/vcf.h>

#include "Variant.hh"

namespace variant
{
    /** store regions for quantification in a file */
    class QuantifyRegions
    {
    public:
        QuantifyRegions();
        ~QuantifyRegions();
        /** load named regions
         *
         *  Each name must give a region name and a bed file:
         *
         *  FP:fp.bed
         *
         */
        void load(std::vector<std::string> const & rnames, bool fixchr=false);

        /**
         * Returns true if regions were loaded.
         *
         * Note that this will also return true when an empty bed file was loaded.
         * This is intentional to distinguish the case where we don't have confident
         * regions (everything unknown is a FP) from the one where the confident
         * region file is empty (every FP is unknown).
         */
        bool hasRegions() const;

        /** add Regions annotation to a record
         *
         * Records must be passed in sorted order.
         *
         */
        void annotate(bcf_hdr_t * hdr, bcf1_t * record);
    private:
        struct QuantifyRegionsImpl;
        std::unique_ptr<QuantifyRegionsImpl> _impl;
    };
}

#endif //HAPLOTYPES_QUANTIFYREGIONS_HH_HH
