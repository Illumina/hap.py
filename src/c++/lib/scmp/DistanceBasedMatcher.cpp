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
 *  \brief Match alleles via normalisation + hashing
 *
 * \file DistanceBasedMatcher.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "DistanceBasedMatcher.hh"
#include "HapSetMatcherImpl.hh"
#include "helpers/IntervalBuffer.hh"

namespace variant
{

    DistanceBasedMatcher::DistanceBasedMatcher(std::string const & reference,
                                 std::string const & chr,
                                 int max_distance_) :
        HapSetMatcher(reference, chr, 1, 0),
        max_distance(max_distance_)
    { }

    DistanceBasedMatcher::~DistanceBasedMatcher() {}

    /**
     * update / optimize assignments
     */
    size_t DistanceBasedMatcher::optimize(HapAssignment & assignment)
    {
        _impl->_updateSortedVariants();
        reset(assignment);

        intervals::IntervalBuffer ivb;
        for(auto const & sv : _impl->sorted_variants)
        {
            const auto v_id = sv.second;
            auto const & v = _impl->variants[v_id];
            const auto v_start = std::min(v.data.start, v.data.end);
            const auto v_end = std::max(v.data.start, v.data.end);
            const int v_lane = (v.side == HapSetMatcherImpl::Variant::LEFT) ? 0 : 1;
            ivb.addInterval(v_start, v_end, v_lane);
        }

        for(auto const & sv : _impl->sorted_variants)
        {
            const auto v_id = sv.second;
            auto const & v = _impl->variants[v_id];
            const auto v_start = std::min(v.data.start, v.data.end) - max_distance;
            const auto v_end = std::max(v.data.start, v.data.end) + max_distance;
            const int v_lane_to_check = (v.side == HapSetMatcherImpl::Variant::LEFT) ? 1 : 0;
            if(ivb.hasOverlap(v_start, v_end, v_lane_to_check))
            {
                assignment.variant_assignments[v_id] = 1;
            }
            else
            {
                assignment.variant_assignments[v_id] = 0;
            }
        }
        return _impl->scoreAssignment(assignment);
    }

    /**
     * update / optimize assignments
     *
     * return < true/false if assignment produces match, number of matched variants >
     * Total number of variants is in assignment.variant_assignments.size()
     */
    std::pair<bool, size_t> DistanceBasedMatcher::checkAndScore(HapAssignment const & assignment)
    {
        const bool valid = true;
        const size_t score = _impl->scoreAssignment(assignment);
        return std::make_pair(valid, score);
    }
}
