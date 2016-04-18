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
 * Variant counting wrapper implementation
 *
 * \file BlockQuantifyImpl.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_BLOCKQUANTIFYIMPL_HH
#define HAPLOTYPES_BLOCKQUANTIFYIMPL_HH

#include "helpers/StringUtil.hh"
#include "helpers/BCFHelpers.hh"
#include "Variant.hh"
#include "Fasta.hh"
#include "Error.hh"

#include <htslib/hts.h>

#include <array>
#include <list>
#include <htslib/vcf.h>
#include <thread>

#include "BlockQuantify.hh"
#include "helpers/Roc.hh"

namespace variant
{
    struct BlockQuantify::BlockQuantifyImpl
    {
        ~BlockQuantifyImpl();

        bcf_hdr_t * hdr;
        FastaFile const & ref_fasta;
        std::unique_ptr<FastaFile> fasta_to_use;

        typedef std::list<std::string> samplenames_t;
        typedef std::map<std::string, VariantStatistics> count_map_t;
        typedef std::list<bcf1_t *> variantlist_t;
        typedef std::map<std::string, roc::Roc> rocmap_t;
        typedef std::set<std::string> filterset_t;

        count_map_t count_map;
        variantlist_t variants;
        samplenames_t samples;
        rocmap_t rocs;
        filterset_t filters_to_ignore;

        bool count_unk;
        bool output_vtc;
        bool count_homref;
        bool extended_counts;
    };
}

#endif //HAPLOTYPES_BLOCKQUANTIFYIMPL_HH
