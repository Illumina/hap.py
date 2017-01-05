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
 * \file HapSetMatcher.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "HapSetMatcher.hh"
#include "HapSetMatcherImpl.hh"
#include "Haplotype.hh"
#include "helpers/Popcount.hh"

#include <map>
#include <unordered_map>

namespace variant {

    HapSetMatcher::HapSetMatcher(std::string const & ref,
                                 std::string const & chr,
                                 int n_haps,
                                 int n_enum
    ) : _impl(new HapSetMatcherImpl())
    {
        _impl->ref = ref;
        _impl->chr = chr;
        assert(n_haps > 0 && n_haps <= 8);
        _impl->n_haps = n_haps;
        _impl->n_enum = n_enum;
    }

    HapSetMatcher::~HapSetMatcher() {}

    size_t HapSetMatcher::addLeft(RefVar const & var, int8_t copies, int64_t ps)
    {
        HapSetMatcherImpl::Variant av;
        av.data = var;
        av.side = HapSetMatcherImpl::Variant::LEFT;
        if(ps < 0)
        {
            av.setCopies(copies, _impl->n_haps);
        }
        else
        {
            // one acceptable assignment which gets permuted through PS updates
            av.acceptable_assignments.insert(copies);
        }
        av.ps = ps;
        _impl->variants.push_back(av);
        return _impl->variants.size() - 1;
    }

    size_t HapSetMatcher::addRight(RefVar const & var, int8_t copies, int64_t ps)
    {
        HapSetMatcherImpl::Variant av;
        av.data = var;
        av.side = HapSetMatcherImpl::Variant::RIGHT;
        if(ps < 0)
        {
            av.setCopies(copies, _impl->n_haps);
        }
        else
        {
            // one acceptable assignment which gets permuted through PS updates
            av.acceptable_assignments.insert(copies);
        }
        av.ps = ps;
        _impl->variants.push_back(av);
        return _impl->variants.size() - 1;
    }


    /**
     * update / optimise assignments
     */
    size_t HapSetMatcher::optimize(HapAssignment & assignment)
    {
        if(assignment.variant_assignments.size() != _impl->variants.size())
        {
            reset(assignment);
        }

        int n_enum = _impl->n_enum;

        HapAssignment best_assignment = assignment;
        int64_t best_score = -1;

        struct MaxEnumException {};

        try
        {
            bool ps_left = true;
            while(ps_left)
            {
                bool vars_left = true;
                while(vars_left)
                {
                    if(_impl->isValidAssignment(assignment))
                    {
                        size_t current_score = _impl->scoreAssignment(assignment);
                        if((signed)current_score > best_score)
                        {
                            best_assignment = assignment;
                            best_score = (signed)current_score;
                        }
                    }
                    vars_left = _impl->nextVarAssignment(assignment);
                    if(n_enum-- <= 0)
                    {
                        throw MaxEnumException();
                    }
                }
                ps_left = _impl->nextPSAssignment(assignment);
            }
        }
        catch (MaxEnumException const & )
        {
            // enumeration terminated early?
            best_score = 0;
        }

        assignment = best_assignment;
        return (size_t)best_score;
    }

    /**
     * @return the number of possible assignments
     */
    uint64_t HapSetMatcher::numberOfPossibleAssignments() const
    {
        uint64_t current_count = 1;
        uint64_t ps_count = 0;
        for(auto const & v : _impl->variants)
        {
            uint64_t this_count = v.acceptable_assignments.size();
            if(!v.acceptable_assignments.count(0))
            {
                ++this_count;
            }

            uint64_t next_count = current_count * this_count;
            if(next_count < current_count)
            {
                return std::numeric_limits<uint64_t>::max();
            }
            current_count = next_count;
            if(v.ps >= 0)
            {
                ++ps_count;
            }
        }

        while(ps_count--)
        {
            uint64_t next_count = current_count << 1;
            if(next_count < current_count)
            {
                return std::numeric_limits<uint64_t>::max();
            }
            current_count = next_count;
        }
        return current_count;
    }


    /**
     * reset/initialize an assignment
     */
    void HapSetMatcher::reset(HapAssignment & assignment)
    {
        _impl->resetAssignment(assignment);
    }

    /**
     * update / optimize assignments
     *
     * return < true/false if assignment produces match, number of matched variants >
     * Total number of variants is in assignment.variant_assignments.size()
     */
    std::pair<bool, size_t> HapSetMatcher::checkAndScore(HapAssignment const & assignment)
    {
        const bool valid = _impl->isValidAssignment(assignment);
        const size_t score = valid ? _impl->scoreAssignment(assignment) : 0;
        return std::make_pair(valid, score);
    }
}
