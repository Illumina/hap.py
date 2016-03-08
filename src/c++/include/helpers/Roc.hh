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
 * \file Roc.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */
#ifndef HAPLOTYPES_ROC_HH
#define HAPLOTYPES_ROC_HH

#include <memory>
#include <vector>

namespace roc
{
    template <typename E>
    static constexpr typename std::underlying_type<E>::type to_underlying(E e) {
        return static_cast<typename std::underlying_type<E>::type>(e);
    }
    enum class DecisionType : int { FN, TP, FP, UNK, N, SIZE };
    const int NDecisionTypes = to_underlying(DecisionType::SIZE);
    extern const char * DecisionTypes[NDecisionTypes];

    struct Observation
    {
        double level;
        DecisionType dt;
        uint64_t n;
    };

    struct Level
    {
        Level() : level(-1) { for(int x = 0; x < to_underlying(DecisionType::SIZE); ++x) { counts[x] = 0; } }
        double level;
        uint64_t counts[to_underlying(DecisionType::SIZE)];

        /** add single observation */
        void addObs(Observation const & obs)
        {
            counts[to_underlying(obs.dt)] += obs.n;
        }
    };

    class Roc {
    public:
        Roc();
        ~Roc();
        Roc(Roc const & rhs);
        Roc(Roc && rhs);
        Roc & operator=(Roc && rhs);

        Roc & operator=(Roc const& rhs) = delete;

        // add observations from second ROC
        void add(Roc const &rhs);
        void add(Observation const &rhs);

        void getLevels(std::vector<Level> & target, bool reversed=false) const;
    private:
        struct RocImpl;
        std::unique_ptr<RocImpl> _impl;
    };
}


#endif //HAPLOTYPES_ROC_HH
