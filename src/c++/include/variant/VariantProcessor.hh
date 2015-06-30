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
 *
 * \file VariantProcessor.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"

namespace variant
{

/** Process variants through set of steps */
class VariantProcessor
{
public:
    VariantProcessor();
    VariantProcessor(VariantProcessor const & rhs);
    ~VariantProcessor();
    VariantProcessor const & operator=(VariantProcessor const & rhs);

    /** set up processing */
    void addStep(AbstractVariantProcessingStep &, bool prepend=false);

    /** process a Variant Reader */
    void setReader(VariantReader & input, VariantBufferMode mode, int64_t param=0);

    /**
     * @brief Rewind / set region to read
     *
     * @param chr chromosome/contig name
     * @param startpos start position on chr (or -1 for entire chr)
     *
     *
     * Note that this involves clearing internal buffers, which might be slow if
     * the lots of variants are in the pipeline still.
     */
    void rewind(const char * chr=NULL, int64_t startpos=-1);

    /** Variant output **/
    /**
     * @brief Return variant block at current position
     **/
    Variants & current();

    /**
     * @brief Advance one line
     * @return true if a variant was retrieved, false otherwise
     */
    bool advance();

    /** put back a variant to the last stage */
    void putBack(Variants const & v);
private:
    struct VariantProcessorImpl;
    VariantProcessorImpl * _impl;
};

}
