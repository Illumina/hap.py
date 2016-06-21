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
    /**
     * Decision types:
     *
     * on a ROC, we change variants depending if their level is above (PASS) or below (FAIL) the current threshold
     *
     * reported     PASS    FAIL    comment
     * -------------------------------------------------------
     * FN           FN      FN      counted against recall
     * TP           TP      FN      counted for recall
     * FN2          FN2     FN2     FNs in query representation
     * TP2          TP2     FN2     TPs in query representation
     * FP           FP      N       counted against precision
     * UNK          UNK     N       counted as NA
     * N            N       N       ignored / filtered
     *
     */
    enum class DecisionType : int { FN, TP, FN2, TP2, FP, UNK, N, SIZE };
    const int NDecisionTypes = to_underlying(DecisionType::SIZE);
    extern const char * DecisionTypes[NDecisionTypes];

    /**
     * Flags for observations.
     */

    /** match / mismatch types */
    static const uint64_t OBS_FLAG_LM = 0x0001;
    static const uint64_t OBS_FLAG_AM = 0x0002;
    static const uint64_t OBS_FLAG_GM = 0x0004;

    /** het / hom / hetalt */
    static const uint64_t OBS_FLAG_HET = 0x0020;
    static const uint64_t OBS_FLAG_HOMALT = 0x0040;
    static const uint64_t OBS_FLAG_HETALT = 0x0080;

    /** various BI flags */

    /** SNP TI / TV */
    static const uint64_t OBS_FLAG_TI = 0x0008;
    static const uint64_t OBS_FLAG_TV = 0x0010;

    /** indel lengths */
    static const uint64_t OBS_FLAG_I1_5    = 0x0100;
    static const uint64_t OBS_FLAG_I6_15   = 0x0200;
    static const uint64_t OBS_FLAG_I16_PLUS = 0x0400;
    static const uint64_t OBS_FLAG_D1_5    = 0x0800;
    static const uint64_t OBS_FLAG_D6_15   = 0x1000;
    static const uint64_t OBS_FLAG_D16_PLUS = 0x2000;
    static const uint64_t OBS_FLAG_C1_5    = 0x4000;
    static const uint64_t OBS_FLAG_C6_15   = 0x8000;
    static const uint64_t OBS_FLAG_C16_PLUS = 0x10000;

    uint64_t makeObservationFlags(const std::string & bk, const std::string & bi, const std::string & lt);

    struct Observation
    {
        double level;
        DecisionType dt;
        uint64_t n;
        uint64_t flags;
    };

    struct Level
    {
        Level() : level(-1), _fp_gt(0), _fp_al(0)
        {
            for(int x = 0; x < to_underlying(DecisionType::SIZE); ++x) { counts[x] = 0; }
        }

        double level;
        uint64_t counts[to_underlying(DecisionType::SIZE)];

        /** these counts are added from observations via flags */
        uint64_t _fp_gt;
        uint64_t _fp_al;

        /** add single observation */
        void addObs(Observation const & obs)
        {
            counts[to_underlying(obs.dt)] += obs.n;
            if(obs.dt == DecisionType::FP)
            {
                // count "special" mismatches
                if(obs.flags & OBS_FLAG_AM)
                {
                    _fp_gt += obs.n;
                }
                if(obs.flags & OBS_FLAG_LM)
                {
                    _fp_al += obs.n;
                }
            }
        }

        uint64_t tp() const { return counts[to_underlying(DecisionType::TP)]; }
        uint64_t tp2() const { return counts[to_underlying(DecisionType::TP2)]; }
        uint64_t fp() const { return counts[to_underlying(DecisionType::FP)]; }
        uint64_t fp_gt() const { return _fp_gt; }
        uint64_t fp_al() const { return _fp_al; }
        uint64_t fn() const { return counts[to_underlying(DecisionType::FN)]; }
        uint64_t fn2() const { return counts[to_underlying(DecisionType::FN2)]; }
        uint64_t unk() const { return counts[to_underlying(DecisionType::UNK)]; }
        uint64_t n() const { return counts[to_underlying(DecisionType::N)]; }

        uint64_t tp(uint64_t val) { counts[to_underlying(DecisionType::TP)] = val; return val; }
        uint64_t tp2(uint64_t val) { counts[to_underlying(DecisionType::TP2)] = val; return val; }
        uint64_t fp(uint64_t val) { counts[to_underlying(DecisionType::FP)] = val; return val; }
        uint64_t fp_gt(uint64_t val) { _fp_gt = val; return val; }
        uint64_t fp_al(uint64_t val) { _fp_al = val; return val; }
        uint64_t fn(uint64_t val) { counts[to_underlying(DecisionType::FN)] = val; return val; }
        uint64_t fn2(uint64_t val) { counts[to_underlying(DecisionType::FN2)] = val; return val; }
        uint64_t unk(uint64_t val) { counts[to_underlying(DecisionType::UNK)] = val; return val; }
        uint64_t n(uint64_t val) { counts[to_underlying(DecisionType::N)] = val; return val; }

        uint64_t totalTruth() const { return tp() + fn(); }
        // total variants reported by query at this level.
        uint64_t totalQuery() const { return tp2() + fp() + unk(); }

        /** recall relative to truth */
        double recall() const
        {
            const uint64_t total_truth = totalTruth();
            return (total_truth == 0) ? 0 : static_cast<double>(tp()) / (total_truth);
        }

        /** relative precision of query */
        double precision() const
        {
            const uint64_t total_query = totalQuery();
            return (total_query == 0) ? 0 : static_cast<double>(tp2()) / (tp2() + fp());
        }

        /** fraction of query calls that are UNK */
        double na() const
        {
            const uint64_t total_query = totalQuery();
            return (total_query == 0) ? 0 : static_cast<double>(unk()) / (total_query);
        }

        /** F1 score: https://en.wikipedia.org/wiki/F1_score */
        double fScore() const
        {
            const double p = precision();
            const double r = recall();
            return 2.0 * p*r / (p + r);
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

        void getLevels(std::vector<Level> & target, double roc_delta=0, uint64_t flag_mask=0) const;
        Level getTotals(uint64_t flag_mask=0) const;
    private:
        struct RocImpl;
        std::unique_ptr<RocImpl> _impl;
    };
}


#endif //HAPLOTYPES_ROC_HH
