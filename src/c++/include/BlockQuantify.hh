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
 * Count variants
 *
 * \file BlockQuantify.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_BLOCKQUANTIFY_HH
#define HAPLOTYPES_BLOCKQUANTIFY_HH

#include <memory>
#include <map>
#include <string>

#include <htslib/vcf.h>

#include "Variant.hh"
#include "IntervalTree.h"

namespace variant {

    /** store regions for quantification in a file */
    struct QuantifyRegions
    {
        typedef std::map<std::string, IntervalTree<std::string, int64_t>> regionmap_t;
        void load(std::vector<std::string> const & rnames);
        regionmap_t regions;
    };

    class BlockQuantify {
    public:
        // hdr must be destroyed externally, the reason it's not const is
        // that htslib doesn't do const
        explicit BlockQuantify(bcf_hdr_t * hdr,
                               QuantifyRegions const & qregions,
                               std::string const & ref_fasta,
                               bool output_vtc,
                               bool count_homref);
        ~BlockQuantify();

        BlockQuantify(BlockQuantify && rhs);
        BlockQuantify & operator=(BlockQuantify && rhs);

        BlockQuantify(BlockQuantify const& rhs) = delete;
        BlockQuantify & operator=(BlockQuantify const& rhs) = delete;

        // add a bcf record (will take ownership of this record)
        void add(bcf1_t * v);

        // count saved variants
        void count();

        // result output
        std::map<std::string, VariantStatistics> const & getCounts() const;
        std::list<bcf1_t*> const & getVariants();
    protected:
        void count_variants(bcf1_t * v);

        /** save lookups when reading from the same chromosome */
        void start_chr(std::string const & chr);
    private:
        struct BlockQuantifyImpl;
        std::unique_ptr<BlockQuantifyImpl> _impl;
    };

}

#endif //HAPLOTYPES_BLOCKQUANTIFY_HH
