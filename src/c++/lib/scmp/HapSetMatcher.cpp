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
#include "Haplotype.hh"
#include "helpers/Popcount.hh"

#include <map>
#include <unordered_map>

namespace variant {

    static const int MAX_ACCEPTABLE_ASSIGNMENTS = 2;

    struct HapSetMatcher::HapSetMatcherImpl
    {
        HapSetMatcherImpl(std::string const & ref) : ref(ref.c_str()) {}
        struct Variant
        {
            RefVar data;
            enum {LEFT=0, RIGHT=1} side;
            int64_t ps = -1;
            std::set<int8_t> acceptable_assignments;

            /** initialize acceptable assignments based on number of copies and haps */
            void setCopies(int8_t copies, int n_haps) {
                uint64_t current = 0;
                uint64_t cm = 1;
                assert(copies <= n_haps);
                for(int c = 0; c < copies; ++c)
                {
                    current |= cm;
                    cm <<= 1;
                }
                if(copies < n_haps)
                {
                    cm <<= n_haps - copies;
                }
                while(current < cm)
                {
                    acceptable_assignments.insert((int8_t)current);
                    current = next_permutation(current);
                }
            }

            static std::vector< std::set<int8_t> > precomputed_assignments;
        };


        void resetAssignment(HapAssignment & a)
        {
            a.ps_assignments[0].clear();
            a.ps_assignments[1].clear();
            a.variant_assignments.resize(variants.size());

            uint64_t default_ps_assignment = 0;
            uint64_t current_assignment = 1;
            for(int hap = 0; hap < n_haps; ++hap)
            {
                default_ps_assignment |= current_assignment;
                current_assignment <<= 9; // 8 -> next slot, 1 => move mask
            }

            for(size_t vid = 0; vid < variants.size(); ++vid)
            {
                auto const & v = variants[vid];
                if(v.ps >= 0)
                {
                    a.ps_assignments[v.side][v.ps] = default_ps_assignment;
                }
                if(v.acceptable_assignments.empty())
                {
                    a.variant_assignments[vid] = 0;
                }
                else
                {
                    a.variant_assignments[vid] = *(v.acceptable_assignments.cbegin());
                }
            }
        }

        /**
         * Test if an assignment produces matching haplotypes
         * @param a assignment of variants to haplotypes
         * @return match status -- true when all haplotype pairs match
         */
        bool isValidAssignment(HapAssignment const & a)
        {
            std::vector<haplotypes::Haplotype> haps;
            haplotypes::Haplotype hap(chr.c_str(), ref.getFilename().c_str());
            haps.resize((unsigned long) (2 * n_haps), hap);

            if(sorted_variants.size() != variants.size())
            {
                sorted_variants.clear();
                size_t vid = 0;
                for(auto const & v : variants)
                {
                    sorted_variants.insert(std::make_pair(v.data.start, vid));
                    ++vid;
                }
            }

            int64_t min_pos = std::numeric_limits<int64_t>::max();
            int64_t max_pos = 0;

            for(auto const & v_id : sorted_variants)
            {
                const auto & v = variants[v_id.second];

                min_pos = std::min(min_pos, v.data.start);
                min_pos = std::min(min_pos, v.data.end);
                max_pos = std::max(max_pos, v.data.start);
                max_pos = std::max(max_pos, v.data.end);

                assert(a.variant_assignments.size() > v_id.second);
                int8_t assignment = a.variant_assignments[v_id.second];
                if (v.ps >= 0)
                {
                    auto const &ps_assignment = a.ps_assignments[(int) v.side];
                    auto ps_it = ps_assignment.find(v.ps);
                    if (ps_it == ps_assignment.end())
                    {
                        // unassigned PS
                        return false;
                    }

                    // reassign according to ps
                    int64_t new_assignments = ps_it->second;
                    int8_t updated_assignment = 0;
                    while(assignment > 0)
                    {
                        updated_assignment |= (int8_t )new_assignments;
                        new_assignments >>= 8;
                        assignment >>= 1;
                    }
                    assignment = updated_assignment;
                }

                int assigned_hap = 0;
                while(assignment && assigned_hap < n_haps)
                {
                    if(assignment & 1)
                    {
                        const size_t hpos = (size_t) ((assigned_hap << 1) | v.side);
                        try
                        {
                            haps[hpos].addVar(v.data.start, v.data.end, v.data.alt);
                        }
                        catch(std::runtime_error const & )
                        {
                            // out of order? => this assignment doesn't work.
                            return false;
                        }
                    }
                    assignment >>= 1;
                    ++assigned_hap;
                }
            }

            min_pos = std::max(min_pos - 1, (int64_t)0);
            max_pos = std::max(max_pos + 1, (int64_t)0);

            if(!variants.empty() && min_pos > max_pos)
            {
                return false;
            }

            for(int h = 0; h < n_haps; ++h)
            {
                const std::string h1 = haps[h << 1].seq(min_pos, max_pos);
                const std::string h2 = haps[(h << 1) + 1].seq(min_pos, max_pos);
                if(h1 != h2)
                {
                    return false;
                }
            }

            return true;
        }

        size_t scoreAssignment(HapAssignment const & a)
        {
            size_t score = 0;
            size_t vid = 0;
            for(auto const & v : variants)
            {
                assert(a.variant_assignments.size() > vid);
                auto assignment = a.variant_assignments[vid];
                if(v.acceptable_assignments.count(assignment))
                {
                    score++;
                }
                ++vid;
            }
            return score;
        }

        std::multimap<int64_t, size_t> sorted_variants;
        std::vector<Variant> variants;
        FastaFile ref;
        std::string chr;
        int n_haps = 2;
        int n_enum = 1000;
    };

    HapSetMatcher::HapSetMatcher(std::string const & ref,
                                 std::string const & chr,
                                 int n_haps,
                                 int n_enum
    ) : _impl(new HapSetMatcherImpl(ref))
    {
        _impl->chr = chr;
        assert(n_haps <= 8);
        _impl->n_haps = n_haps;
        _impl->n_enum = n_enum;
    }

    HapSetMatcher::~HapSetMatcher() {}

    size_t HapSetMatcher::addLeft(RefVar const & var, int8_t copies, int64_t ps)
    {
        HapSetMatcherImpl::Variant av;
        av.data = var;
        av.side = HapSetMatcherImpl::Variant::LEFT;
        av.setCopies(copies, _impl->n_haps);
        av.ps = ps;
        _impl->variants.push_back(av);
        return _impl->sorted_variants.size() - 1;
    }

    size_t HapSetMatcher::addRight(RefVar const & var, int8_t copies, int64_t ps)
    {
        HapSetMatcherImpl::Variant av;
        av.data = var;
        av.side = HapSetMatcherImpl::Variant::RIGHT;
        av.setCopies(copies, _impl->n_haps);
        av.ps = ps;
        _impl->variants.push_back(av);
        return _impl->sorted_variants.size() - 1;
    }


    /**
     * update / optimise assignments
     */
    void HapSetMatcher::optimize(HapAssignment & assignment)
    {

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
