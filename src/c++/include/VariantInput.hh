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
 * \brief Wrapper class that creates a Variant Processor and various processing steps
 *
 * \file VariantInput.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"
#include "variant/VariantProcessor.hh"

#include <list>

namespace variant {
    class VariantInput {
    public:
        VariantInput(
            const char * _ref_fasta,
            bool leftshift = false,
            bool refpadding = true,
            bool trimalleles = false,
            bool splitalleles = false,
            int mergebylocation = 0,
            bool uniqalleles = false,
            bool calls_only = true,
            bool homref_split = false,
            bool primitive_split = false,
            bool homref_output = false,
            int64_t leftshift_limit = -1,
            bool collect_raw = false
        );
        ~VariantInput();

        /** Read variants from a block into a list */
        void get(const char * chr, int64_t start, int64_t end, std::list<Variants> & output);

        enum processor_id { variants, homref, raw };
        /** Direct access to variant processor */
        VariantProcessor & getProcessor(processor_id id = variants);

    private:
        // prevent copy construction and assignment
        VariantInput(VariantInput const & ) {}
        VariantInput & operator=(VariantInput const & ) {return *this;}

        struct VariantInputImpl;
        VariantInputImpl * _impl;
    };
}
