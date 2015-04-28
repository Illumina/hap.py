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
 * Aggregate variant calls with alleles starting at the same location
 *
 * \file VariantLocationMap.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"

namespace variant {

typedef std::map<int64_t, variant::Variants> VariantLocationMap;

/**
 * @brief Combine Alleles by Location
 */
VariantLocationMap::iterator addToLocationMap(
    VariantLocationMap & locmap, int sample, 
    RefVar const & rv, bool is_het);

/**
 * @brief Combine alleles in a single Variants record
 */
void addAlleleToVariant(Variants & target, int sample, RefVar const & rv, bool is_het, int ad=-1);

/**
 * @brief Add ref block to variant record
 */
void addRefAlleleToVariant(Variants & target, int sample, int64_t start, int64_t end, bool is_het, int ad=-1);

} // namespace variant
