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
 * \file AlleleMatcher.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "AlleleMatcher.hh"
#include "HapSetMatcherImpl.hh"

namespace variant
{

    AlleleMatcher::AlleleMatcher(std::string const & reference,
                                 std::string const & chr,
                                 int n_enum) :
        HapSetMatcher(reference, chr, 1, n_enum)
    { }

    AlleleMatcher::~AlleleMatcher() {}

    /**
     * update / optimize assignments
     */
    size_t AlleleMatcher::optimize(HapAssignment & assignment)
    {
        _impl->_updateSortedVariants();

        std::multimap<std::string, size_t > left_mapped;

        reset(assignment);
        FastaFile ref(_impl->ref.c_str());

        for(size_t v_id = 0; v_id < _impl->variants.size(); ++v_id)
        {
            auto const & v = _impl->variants[v_id];
            if(v.side == HapSetMatcherImpl::Variant::LEFT)
            {
                RefVar rv = v.data;
                // optional: left-shift
                leftShift(ref, _impl->chr.c_str(), rv);
                trimLeft(ref, _impl->chr.c_str(), rv, false);
                trimRight(ref, _impl->chr.c_str(), rv, false);
                left_mapped.insert(std::make_pair(rv.repr(), v_id));
            }
        }

        for(size_t v_id = 0; v_id < _impl->variants.size(); ++v_id)
        {
            auto const & v = _impl->variants[v_id];
            if(v.side == HapSetMatcherImpl::Variant::RIGHT)
            {
                RefVar rv = v.data;
                leftShift(ref, _impl->chr.c_str(), rv);
                trimLeft(ref, _impl->chr.c_str(), rv, false);
                trimRight(ref, _impl->chr.c_str(), rv, false);
                auto it = left_mapped.find(rv.repr());
                if(it != left_mapped.end())
                {
                    assignment.variant_assignments[it->second] = 1;
                    assignment.variant_assignments[v_id] = 1;
                    left_mapped.erase(it);
                }
                else
                {
                    assignment.variant_assignments[v_id] = 0;
                }
            }
        }
        for(auto & x : left_mapped)
        {
            assignment.variant_assignments[x.second] = 0;
        }

        bool works = _impl->isValidAssignment(assignment);
        if(works)
        {
            return _impl->scoreAssignment(assignment);
        }
        else
        {
            return 0;
        }
    }

}
