// -*- mode: c++; indent-tabs-mode: nil; -*-
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
 * Store and process observations for ROCs
 *
 * \file Roc.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "helpers/Roc.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <helpers/StringUtil.hh>

namespace roc
{
//    enum class DecisionType : int { FN, TP, FN2, TP2, FP, UNK, N, SIZE };
    const char * DecisionTypes[NDecisionTypes] { "TRUTH.FN", "TRUTH.TP", "QUERY.FN", "QUERY.TP", "QUERY.FP", "QUERY.UNK", "N" };


    uint64_t makeObservationFlags(const std::string & bk, const std::string & bi, const std::string & blt)
    {
        uint64_t result = 0;

        if(bk == "lm")
        {
            result |= OBS_FLAG_LM;
        }
        else if(bk == "am")
        {
            result |= OBS_FLAG_AM;
        }
        else if(bk == "gm")
        {
            result |= OBS_FLAG_GM;
        }

        if(blt == "het")
        {
            result |= OBS_FLAG_HET;
        }
        else if(blt == "hetalt")
        {
            result |= OBS_FLAG_HETALT;
        }
        else if(blt == "homalt")
        {
            result |= OBS_FLAG_HOMALT;
        }

        std::vector<std::string> bis;
        stringutil::split(bi, bis, ",");

        for(auto const & b : bis)
        {
            if(b == "ti")
            {
                result |= OBS_FLAG_TI;
            }
            else if(b == "tv")
            {
                result |= OBS_FLAG_TV;
            }
            else if(b == "tv")
            {
                result |= OBS_FLAG_TV;
            }
            else if(b == "i1_5")
            {
                result |= OBS_FLAG_I1_5;
            }
            else if(b == "i6_15")
            {
                result |= OBS_FLAG_I6_15;
            }
            else if(b == "i16_plus")
            {
                result |= OBS_FLAG_I16_PLUS;
            }
            else if(b == "d1_5")
            {
                result |= OBS_FLAG_D1_5;
            }
            else if(b == "d6_15")
            {
                result |= OBS_FLAG_D6_15;
            }
            else if(b == "d16_plus")
            {
                result |= OBS_FLAG_D16_PLUS;
            }
            else if(b == "c1_5")
            {
                result |= OBS_FLAG_C1_5;
            }
            else if(b == "c6_15")
            {
                result |= OBS_FLAG_C6_15;
            }
            else if(b == "c16_plus")
            {
                result |= OBS_FLAG_C16_PLUS;
            }
        }

        return result;
    }


    struct Roc::RocImpl
    {
        std::vector<Observation> obs;
    };

    Roc::Roc() : _impl(new RocImpl()) { }
    Roc::~Roc() { }
    Roc::Roc(Roc const & rhs) : _impl(new RocImpl()) { add(rhs); }
    Roc::Roc(Roc && rhs) : _impl(std::move(rhs._impl)) { }

    Roc & Roc::operator=(Roc && rhs)
    {
        _impl = std::move(rhs._impl);
        return *this;
    }

    // add observations from second ROC
    void Roc::add(Roc const &rhs)
    {
        _impl->obs.insert(_impl->obs.end(), rhs._impl->obs.cbegin(), rhs._impl->obs.cend());
    }

    void Roc::add(Observation const &rhs)
    {
        _impl->obs.push_back(rhs);
    }

    Level Roc::getTotals(uint64_t flag_mask) const
    {
        Level last;
        for(auto const & x : _impl->obs)
        {

            if(flag_mask > 0 && ((x.flags & flag_mask) != flag_mask))
            {
                // skip if required flags are not set
                continue;
            }
            last.addObs(x);
        }
        last.level = std::numeric_limits<double>::quiet_NaN();
        return last;
    }

    void Roc::getLevels(std::vector<Level> & output, double roc_delta, uint64_t flag_mask) const
    {
        std::sort(_impl->obs.begin(), _impl->obs.end(),
                  [](Observation const & o1, Observation const & o2) -> bool {
                      return o1.level < o2.level;
                  });

        Level last;
        std::vector<Level> target;
        last.level = std::numeric_limits<double>::quiet_NaN();
        for(auto const & x : _impl->obs)
        {
            if(flag_mask > 0 && ((x.flags & flag_mask) != flag_mask))
            {
                // skip if required flags are not set
                continue;
            }
            last.addObs(x);
            last.level = x.level;
            // this is probably not the best way to do this, but it works
            // because we print things to a text file in the end
            target.push_back(last);
        }

        for(auto & x : target)
        {
            // TPs above or on level
            const uint64_t tp = last.tp() - x.tp();
            // FPs above or on level
            const uint64_t fp = last.fp() - x.fp();
            const uint64_t fp_gt = last.fp_gt() - x.fp_gt();
            const uint64_t fp_al = last.fp_al() - x.fp_al();
            // UNKs above or on level
            const uint64_t unk = last.unk() - x.unk();

            // FN = last level FNs + TPs below level
            const uint64_t fn = last.fn() + x.tp();

            // TPs above level
            const uint64_t tp2 = last.tp2() - x.tp2();
            // FN = FN at or below this level + TPs below level
            const uint64_t fn2 = last.fn2() + x.tp2();

            // N = n + fp below level + unk below level
            const uint64_t n = last.n() + x.fp() + x.unk();

            x.fn(fn);
            x.fn2(fn2);
            x.n(n);
            x.tp(tp);
            x.tp2(tp2);
            x.fp(fp);
            x.fp_gt(fp_gt);
            x.fp_al(fp_al);
            x.unk(unk);
        }

        if(target.empty())
        {
            return;
        }
        double previous_level = target.front().level;
        bool is_first = true;
        bool push_now = false;
        for(auto const & x : target)
        {
            if(roc_delta < std::numeric_limits<double>::epsilon())
            {
                const std::string cur = std::to_string(x.level);
                const std::string prev = std::to_string(previous_level);
                push_now = cur != prev;
            }
            else
            {
                if(is_first)
                {
                    push_now = true;
                    is_first = false;
                }
                else
                {
                    push_now = fabs(x.level - previous_level) > roc_delta;
                }
            }
            if(push_now)
            {
                output.push_back(x);
                previous_level = x.level;
            }
        }
    }
}
