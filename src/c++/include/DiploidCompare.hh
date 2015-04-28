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
 * \brief Diploid partially phased haplotype comparison
 *  
 * \file DiploidCompare.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "GraphReference.hh"
#include "DiploidReference.hh"
#include "RefVar.hh"
#include "HaploCompare.hh"

#include "DiploidComparisonResult.hh"

namespace haplotypes
{

/**
 * @brief Diploid Haplotype comparison
 * @details Compare Variants assuming 
 * 
 */
struct DiploidCompareImpl;
class DiploidCompare
{
public:
    DiploidCompare(const char * ref_fasta);
    DiploidCompare(DiploidCompare const & );
    ~DiploidCompare();
    DiploidCompare const & operator=(DiploidCompare const & );

    /**
     * Set the maximum number of haplotype blocks to enumerate before aborting
     */
    void setMaxHapEnum(int nhap=4096);

    /**
     * Enable / disable the alignment step to find the best approximately matching
     * haplotypes (i.e. stop after match/mismatch status have been established)
     */
    void setDoAlignments(bool doAlignments = true);
    bool getDoAlignments();

    /**
     * Enable / disable the multiple alignment step to find canonical representations
     * of shared and unique variation.
     */
    void setDoSharedVar(bool doSharedVar = false);
    bool getDoSharedVar();

    /**
     * @brief Set the region to compare in and reset the enumeration.
     * 
     * Read variants from list rather than the graph references
     */
    void setRegion(const char * chr, int64_t start, int64_t end,
                   std::list<variant::Variants> const & vars, int ix1, int ix2);

    /**
     * @brief return comparison outcome after setRegion
     */
    DiploidComparisonResult const & getResult();

private:
    DiploidCompareImpl * _impl;
};


} // namespace haplotypes

