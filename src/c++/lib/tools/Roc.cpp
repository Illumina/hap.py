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

namespace roc
{
    const char * DecisionTypes[NDecisionTypes] { "FN", "TP", "FP", "UNK", "N" };

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

    void Roc::getLevels(std::vector<Level> & target, bool reversed) const
    {
        if(!reversed)
        {
            std::sort(_impl->obs.begin(), _impl->obs.end(),
                      [](Observation const & o1, Observation const & o2) -> bool {
                          return o1.level < o2.level;
                      });
        }
        else
        {
            std::sort(_impl->obs.begin(), _impl->obs.end(),
                      [](Observation const & o1, Observation const & o2) -> bool {
                          return o2.level < o1.level;
                      });
        }


        Level last;
        last.level = std::numeric_limits<double>::quiet_NaN();
        for(auto const & x : _impl->obs)
        {
            last.addObs(x);
            // this is probably not the best way to do this, but it works
            // because we print things to a text file in the end
            if(!std::isnan(last.level) && std::to_string(last.level) != std::to_string(x.level))
            {
                target.push_back(last);
            }
            last.level = x.level;
        }

        for(auto & x : target)
        {
            // TPs above level
            const uint64_t tp = last.tp() - x.tp();
            // FPs above level
            const uint64_t fp = last.fp() - x.fp();
            // UNKs above level
            const uint64_t unk = last.unk() - x.unk();

            // FN = FN + TPs below level
            const uint64_t fn = last.fn() + x.tp();

            // TPs above level
            const uint64_t tp2 = last.tp2() - x.tp2();
            // FN = FN + TPs below level
            const uint64_t fn2 = last.fn2() + x.tp2();

            // FN = FN + FPs below level + UNK below level
            const uint64_t n = last.n() + x.fp() + x.unk();

            x.fn(fn);
            x.fn2(fn2);
            x.n(n);
            x.tp(tp);
            x.tp2(tp2);
            x.fp(fp);
            x.unk(unk);
        }

    }
}
