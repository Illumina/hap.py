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
 * Compares alleles in blocks
 *
 * \file BlockAlleleCompare.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_BLOCKALLELECOMPARE_HH_HH
#define HAPLOTYPES_BLOCKALLELECOMPARE_HH_HH

#include <memory>
#include <map>
#include <string>

#include <htslib/vcf.h>

#include "Variant.hh"
#include "json/json.h"

namespace variant {

    class BlockAlleleCompare {
    public:
        /**
         * hdr must be destroyed externally, the reason it's not const is
         * that htslib doesn't do const mostly.
         *
         * @param hdr a bcf header
         * @param ref_fasta reference fasta file for trimming and counting
         * @param qq_name name of feature to use for QQ. Can be QUAL, an INFO or a FORMAT field
         */
        explicit BlockAlleleCompare(bcfhelpers::p_bcf_hdr hdr,
                                    FastaFile const & ref_fasta,
                                    std::string const & qq_name);
        virtual ~BlockAlleleCompare();

        BlockAlleleCompare(BlockAlleleCompare && rhs);
        BlockAlleleCompare & operator=(BlockAlleleCompare && rhs);

        BlockAlleleCompare(BlockAlleleCompare const& rhs) = delete;
        BlockAlleleCompare & operator=(BlockAlleleCompare const& rhs) = delete;

        /**
         * operators so we can sort blocks
         */
        bool operator==(BlockAlleleCompare const & rhs) const;
        bool operator<(BlockAlleleCompare const & rhs) const;

        enum ComparisonMode {ALLELES, DISTANCE, ENUMERATE_DIPLOID};

        /**
         * Set the comparison mode, and reset the
         * comparison parameters to default.
         */
        void setComparisonMode(ComparisonMode mode);

        /**
         * Set comparison parameters
         * @param params parameters for the comparison method
         */
        void setComparisonParameters(Json::Value const & params);

        /**
         * Add a BCF record. Will duplicate the record and keep the copy
         * @param v bcf record
         */
        void add(bcf1_t * v);

        /**
         * Compare all buffered variants
         */
        virtual void run();

        /**
         * Output all buffered variants into a file
         * @param output HTS file to write to
         */
        void output(htsFile * output);

        /**
         * BlockAlleleCompare::run doesn't write a header such that
         * results can be concatenated. This function updates the header once in the beginning
         * @param hdr header to update
         */
        static void updateHeader(bcf_hdr_t * hdr);
    protected:
        struct BlockAlleleCompareImpl;
        std::unique_ptr<BlockAlleleCompareImpl> _impl;
    };

}


#endif //HAPLOTYPES_BLOCKALLELECOMPARE_HH_HH
