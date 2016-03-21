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
#include "Fasta.hh"

#include <json/json.h>
#include <list>
#include <vector>

#include <htslib/vcf.h>

namespace variant
{
typedef std::vector<std::string> StrVec;

struct VariantStatisticsImpl;
class VariantStatistics
{
public:
    VariantStatistics(const char *, bool _tmp=false) = delete;
    VariantStatistics(std::string, bool _tmp=false) = delete;
    VariantStatistics(FastaFile const & ref_fasta, bool count_homref=false);
    VariantStatistics(VariantStatistics const & rhs);
    ~VariantStatistics();
    VariantStatistics const & operator=(VariantStatistics const & rhs);

    /**
     * @brief Read and write from/to JSON files
     */
    void read(Json::Value const & );
    Json::Value write() const;

    /**
     * @brief add counts
     */
    void add(VariantStatistics const & rhs);
    /** add for variants and alleles can return the types that were added to */
    void add(bcf_hdr_t * hdr, bcf1_t * rhs, int sample, int ** types = NULL, int * ntypes = NULL,
             std::set<std::string> ** extra_counts_seen = NULL);
    void add(Variants const & rhs, int sample, int ** types = NULL, int * ntypes = NULL);
    void add(const char * chr, RefVar const & rhs, int ** types = NULL, int * ntypes = NULL);

    /** resolve types to strings */
    static std::string type2string(int type);

    /** resolve types to strings */
    static int string2type(const char * str);

    void getExtraCountNames(StrVec& extraCountNames);
    size_t extraCount(const std::string& extraCountName);

    /* abbreviate a set of extra counts that were observed for VCF output */
    static std::string extraCountsToBI(std::set<std::string> const & extra_counts_seen);
private:
    VariantStatisticsImpl * _impl;
};

// variant type names
extern const char * VT_NAMES [];
// encode allele types seen in low 4 bits
extern const uint64_t VT_NOCALL;
extern const uint64_t VT_SNP;
extern const uint64_t VT_INS;
extern const uint64_t VT_DEL;
extern const uint64_t VT_REF;


// calltype names -- look up via calltype >> 4
extern const char * CT_NAMES [];
// this is encoded in the high 4 bits 
extern const uint64_t CT_NUCLEOTIDES;     // count nucleotides
extern const uint64_t CT_ALLELES;         // count alleles
extern const uint64_t CT_HOMREF;          // locations with only one allele seen with copy number 2
extern const uint64_t CT_HET;             // locations with one ref and one alt allele seen with copy number 1
extern const uint64_t CT_HETALT;          // locations with two alt alleles seen with copy number 1 each
extern const uint64_t CT_HEMI;            // locations with only one allele seen with copy number 1
extern const uint64_t CT_AMBI;            // locations with more than two alleles
extern const uint64_t CT_HALFCALL;        // locations with one allele, and one missing call
extern const uint64_t CT_NOCALL;          // locations with one allele, and one missing call
extern const uint64_t CT_HOMALT;          // locations with one allele, and one missing call
extern const uint64_t CT_UNKNOWN;         // unknown locations / counts


}
