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
 * VariantProcessor to remove unused alleles
 *
 * \file VariantAlleleRemover.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"

namespace variant
{

/**
 * @brief Remove unused alleles
 */
void trimAlleles(Variants & vars);

class VariantAlleleRemover : public AbstractVariantProcessingStep
{
public:
    VariantAlleleRemover() : firstone(false) {}
    ~VariantAlleleRemover() {}

    /** Variant input **/
    /** enqueue a set of variants */
    void add(Variants const & vs) {
        if (buffer.empty())
        {
            firstone = true;
        }
        buffer.push_back(vs);
        trimAlleles(buffer.back());
    }

    /** Variant output **/
    /**
     * @brief Return variant block at current position
     **/
    Variants & current() { if( buffer.empty() ) { return tmp; } else { return buffer.front(); }  }

    /**
     * @brief Advance one line
     * @return true if a variant was retrieved, false otherwise
     */
    bool advance() {
        if (firstone && !buffer.empty())
        {
            firstone = false;
            return true;
        }
        else
        {
            if(!buffer.empty())
            {
                buffer.pop_front();
            }
            return !buffer.empty();
        }
    }

    /** empty internal buffer */
    void flush() { buffer.clear(); }

private:
    std::list<Variants> buffer;
    Variants tmp;
    bool firstone;
};

}
