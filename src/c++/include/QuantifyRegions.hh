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
#include <unordered_map>
#include <string>

#include <htslib/vcf.h>

#include "Variant.hh"

namespace variant
{
    /** store regions for quantification in a file */
    class QuantifyRegions
    {
    public:
        explicit QuantifyRegions(std::string const & reference);
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
         * Returns true if confident regions were loaded.
         * (i.e. if regions named "CONF" are present)
         */
        bool hasRegions(std::string const & rname) const;

        /**
         * @return all region names
         */
        std::list<std::string> getRegionNames() const;

        /** add Regions annotation to a record
         *
         * Records must be passed in sorted order.
         *
         */
        void annotate(bcf_hdr_t * hdr, bcf1_t * record);

        /**
         * Get total region sizes in NT
         * @param region_name
         * @return  the region size
         */
        size_t getRegionSize(std::string const & region_name) const;

        /**
         * Get level for a region (0 for top-level, 1 for regions specified in a bed file)
         * @param region_name
         * @return level for the given region
         */
        size_t getRegionLevel(std::string const & region_name) const;

        /**
         * Mark one region to intersect with all others. This will
         * @param base base region (e.g. "CONF")
         */
        void setIntersectRegion(std::string const & base);

        /**
         * Get all intersection / extra counts for a given region. By default, will return one
         * pair with <"Subset.Size", size of the region>
         * @param r the region to query counts for
         * @return a list of counts and names
         */
        std::list< std::pair<std::string, size_t> > getRegionExtraCounts(std::string const & r) const;

    private:

        /**
         * Get region intersection sizes in NT
         * @param region_name base region to intersect all others with (e.g. "CONF")
         * @return  mapping of sizes and overlaps by region name
         */
        std::unordered_map<std::string, size_t> getRegionIntersectionSize(std::string const & region_name) const;

        struct QuantifyRegionsImpl;
        std::unique_ptr<QuantifyRegionsImpl> _impl;
    };
}

#endif //HAPLOTYPES_QUANTIFYREGIONS_HH_HH
