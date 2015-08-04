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
 * \brief Count variant numbers from various intermediate data structures
 *
 * \file VariantStatistics.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"

#include <json/json.h>
#include <list>

namespace variant
{

struct VariantStatisticsImpl;
class VariantStatistics
{
public:
    VariantStatistics(const char * ref_fasta, bool count_homref=false);
    VariantStatistics(VariantStatistics const & rhs);
    ~VariantStatistics();
    VariantStatistics const & operator=(VariantStatistics const & rhs);

    /**
     * @brief Read and write from/to JSON files
     */
    void read(Json::Value const & );
    Json::Value write();

    /**
     * @brief add counts
     */
    void add(VariantStatistics const & rhs);
    /** add for variants and alleles can return the types that were added to */
    void add(Variants const & rhs, int sample, int ** types = NULL, int * ntypes = NULL);
    void add(const char * chr, RefVar const & rhs, int ** types = NULL, int * ntypes = NULL);

    /** resolve types to strings */
    std::string type2string(int type);

    /** resolve types to strings */
    int string2type(const char * str);
private:
    VariantStatisticsImpl * _impl;
};

}
