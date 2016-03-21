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
 * Count variants coming out of xcmp, use hap.py output format
 *
 * \file XCmpQuantify.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_XCMPQUANTIFY_HH
#define HAPLOTYPES_XCMPQUANTIFY_HH

#include "BlockQuantify.hh"
#include "BlockQuantifyImpl.hh"

namespace variant {
    class XCMPQuantify : public BlockQuantify {
    public:
        /**
         * hdr must be destroyed externally, the reason it's not const is
         * that htslib doesn't do const mostly.
         *
         * @param hdr a bcf header
         * @param ref_fasta reference fasta file for trimming and counting
         * @param params string that may contain any of the following:
         *               "output_vtc" add a VTC info field which gives details on what was counted
         *               "count_homref" set to true to also count REF matches
         *               "count_unk" everything that has a missing Regions info tag will become UNK rather
         *                           than FP
         */
        explicit XCMPQuantify(bcf_hdr_t * hdr,
                              FastaFile const & ref_fasta,
                              const std::string & params);
        virtual void updateHeader(bcf_hdr_t * hdr);
    protected:
        virtual void countVariants(bcf1_t * v);
    private:
        // field to use for ROCs (will be translated into QQ format field)
        std::string roc_field;
        int roc_hdr_id;
        bool roc_field_is_info;
        bool roc_field_is_qual;
        // clean the INFO fields (only keep the GA4GH-compliant ones)
        bool clean_info;
    };
}

#endif //HAPLOTYPES_XCMPQUANTIFY_HH
