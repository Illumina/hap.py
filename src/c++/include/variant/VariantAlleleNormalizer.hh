// -*- mode: c++; indent-tabs-mode: nil; -*-
//
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
 * Variant processing step to normalize alleles by left-shifting and trimming
 *
 * \file VariantAlleleNormalizer.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"

namespace variant {

    class VariantAlleleNormalizer : public AbstractVariantProcessingStep
    {
    public:
        VariantAlleleNormalizer();
        VariantAlleleNormalizer(VariantAlleleNormalizer const & );
        ~VariantAlleleNormalizer();
        VariantAlleleNormalizer const & operator=(VariantAlleleNormalizer const & );

        /**
         * @brief Reference for shifting
         *
         */
        void setReference(std::string const & fasta);

        /**
         * @brief Enable / disable reference padding
         *
         * RefVars can be trimmed to include / not include reference bases
         */
        bool getEnableRefPadding() const;
        void setEnableRefPadding(bool padding=false);

        /**
         * @brief Limit / left-most coordinate for shifting
         */
        int64_t getLeftshiftLimit() const;
        void setLeftshiftLimit(int64_t limit=-1);

        /** return homref variants */
        bool getEnableHomrefVariants() const;
        void setEnableHomrefVariants(bool homref=true);

        /** enqueue a set of variants */
        void add(Variants const & vs);

        /**
         * @brief Return variant block at current position
         **/
        Variants & current();

        /**
         * @brief Advance one line
         * @return true if a variant was retrieved, false otherwise
         */
        bool advance();

        /** empty internal buffer */
        void flush();

    private:
        struct VariantAlleleNormalizerImpl;
        VariantAlleleNormalizerImpl * _impl;
    };

}
