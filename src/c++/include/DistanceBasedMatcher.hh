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
 *  \brief Match variants in block if they are closer than max distance
 *
 * \file DistanceBasedMatcher.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once


#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>

#include "HapSetMatcher.hh"


namespace variant
{
    /**
     * Exact matcher for sets of variants assigned to haplotypes
     *
     * Exponential time with cutoff
     */
    class DistanceBasedMatcher : public HapSetMatcher
    {
    public:
        DistanceBasedMatcher(std::string const &reference,
                             std::string const &chr,
                             int max_distance=30);

        ~DistanceBasedMatcher();

        /**
         * update / optimize assignments
         */
        virtual size_t optimize(HapAssignment &assignment);

        virtual std::pair<bool, size_t> checkAndScore(HapAssignment const & assignment);

    private:
        int max_distance;
    };
}
