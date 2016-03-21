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
#include "helpers/Roc.hh"

namespace variant {

    class BlockQuantify;

    /** factory method */
    std::unique_ptr<BlockQuantify> makeQuantifier(
        bcf_hdr_t * hdr,
        FastaFile const & ref_fasta,
        std::string const & type,
        std::string const & params);

    /** list all quantification methods */
    std::list<std::string> listQuantificationMethods();

    class BlockQuantify {
    public:
        /**
         * hdr must be destroyed externally, the reason it's not const is
         * that htslib doesn't do const mostly.
         *
         * @param hdr a bcf header
         * @param ref_fasta reference fasta file for trimming and counting
         */
        explicit BlockQuantify(bcf_hdr_t * hdr,
                               FastaFile const & ref_fasta,
                               std::string const & params);
        virtual ~BlockQuantify();

        BlockQuantify(BlockQuantify && rhs);
        BlockQuantify & operator=(BlockQuantify && rhs);

        BlockQuantify(BlockQuantify const& rhs) = delete;
        BlockQuantify & operator=(BlockQuantify const& rhs) = delete;

        // set up ROC filtering
        // filters can be a comma-separated list of filters to remove, or "*" to remove all
        void rocFiltering(std::string const & filters);

        // add a bcf record (will take ownership of this record)
        void add(bcf1_t * v);

        // count saved variants
        void count();

        // result output
        std::map<std::string, VariantStatistics> const & getCounts() const;
        std::map<std::string, roc::Roc> const & getRocs() const;
        std::list<bcf1_t*> const & getVariants();
        // update a header with new required fields
        virtual void updateHeader(bcf_hdr_t * hdr) = 0;
    protected:
        // overload to implement actual record counting
        virtual void countVariants(bcf1_t * v) = 0;

        // GA4GH-VCF field-based ROC counting
        virtual void rocEvaluate(bcf1_t * v);

        // add ROC decision point
        void addROCValue(std::string const & roc_identifier,
                         roc::DecisionType dt,
                         double level,
                         uint64_t n,
                         bcf1_t * v);

        struct BlockQuantifyImpl;
        std::unique_ptr<BlockQuantifyImpl> _impl;
    };

}

#endif //HAPLOTYPES_BLOCKQUANTIFY_HH
