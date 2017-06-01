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
 *  \brief Optimizer to find optimal allele replay onto haplotypes
 *
 *
 * \file HapSetMatcher.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once


#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>

#include "Fasta.hh"
#include "Variant.hh"


namespace variant
{
    struct HapAssignment
    {
        /**
         * an assignment is a bit-mask, each bit tells us whether the variant
         * is present on a given haplotype
         */

        std::vector<int8_t> variant_assignments;

        // a PS assignment consists of 8 bytes, each re-mapping
        // one of the haplotypes to a new haplotype mask
        // e.g.: 2 haps, ps assignment:
        // 0x 01 02
        //    h2 h1
        //
        // variant assignment 0x01 (het on hap 1) ... becomes with this PS
        // 0x02 (= 1 * 0x02 + 0* 0x01)
        // variant assignment 0x02 (het on hap 2) ... becomes with this PS
        // 0x01 (= 0 * 0x02 + 1* 0x01)
        // variant assignment 0x03 (both on hap1 / hap 2) ... becomes with this PS
        // 0x03 (= 1 * 0x02 + 1* 0x01)
        std::unordered_map<int64_t, uint64_t> ps_assignments[2];
    };

    /**
     * Exact matcher for sets of variants assigned to haplotypes
     *
     * Exponential time with cutoff
     */
    class HapSetMatcher
    {
    public:
        HapSetMatcher(std::string const & reference,
                      std::string const & chr,
                      int n_haps = 2, int n_enum=1 << 20);
        virtual ~HapSetMatcher();

        /**
         * Add variants to left / right side
         * @param var RefVar allele
         * @param copies number of copies / haplotypes this needs to be assigned to if ps == -1,
         *               otherwise the phased assignment
         * @param ps phase set
         * @return identifier for the variant to recover it in the HapAssignment objects
         */
        size_t addLeft(RefVar const & var, int8_t copies, int64_t ps=-1);
        size_t addRight(RefVar const & var, int8_t copies, int64_t ps=-1);

        /**
         * update / optimize assignments
         */
        virtual size_t optimize(HapAssignment & assignment);

        /**
         * @return the number of possible assignments
         */
        uint64_t numberOfPossibleAssignments() const;

        /**
         * update / optimize assignments
         *
         * return < true/false if assignment produces match, number of matched variants >
         * Total number of variants is in assignment.variant_assignments.size()
         */
        virtual std::pair<bool, size_t> checkAndScore(HapAssignment const & assignment);

        /**
         * reset/initialize an assignment
         */
        void reset(HapAssignment & assignment);

    protected:
        struct HapSetMatcherImpl;
        std::unique_ptr<HapSetMatcherImpl> _impl;
    };
}