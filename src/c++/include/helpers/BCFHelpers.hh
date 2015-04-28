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
 *  \brief C++ Helper functions for HTSlib
 *
 *
 * \file BCFHelpers.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "../Variant.hh"

extern "C" {

// GCC warns us about some things in htslib here. We don't care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>

#pragma GCC diagnostic pop
}

namespace bcfhelpers {

/**
 * @brief Set hg19 contig names in header.
 * @param header BCF header
 */
void bcfHeaderHG19(bcf_hdr_t * header);

/**
 * @brief Retrieve an info field as an integer
 * @details [long description]
 * 
 * @param result the default to return if the field is not present
 */
int getInfoInt(bcf_hdr_t * header, bcf1_t * line, const char * field, int result = -1);

/**
 * @brief Read the GT field
 */
void getGT(bcf_hdr_t * header, bcf1_t * line, int isample, int * gt, int & ngt, bool & phased);

/** read GQ(X) -- will use in this order: GQ, GQX, -1 */
void getGQ(const bcf_hdr_t * header, bcf1_t * line, int isample, float & gq);

/** read AD */
void getAD(const bcf_hdr_t * header, bcf1_t * line, int isample, int *ad, int max_ad);

/** read DP(I) -- will use in this order: DP, DPI, -1 */
void getDP(const bcf_hdr_t * header, bcf1_t * line, int isample, int & dp);

} // namespace bcfhelpers
