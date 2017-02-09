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
 * Advanced ROC computations (Ti/Tv, AuC, ...) and ROC output
 *
 * \file RocOutput.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_ROCOUTPUT_HH
#define HAPLOTYPES_ROCOUTPUT_HH

#include "Roc.hh"
#include "QuantifyRegions.hh"

namespace roc
{
    typedef std::map<std::string, roc::Roc> RocMap;

    class ROCOutput
    {
    public:
        /**
         * ROC output class
         *
         * @param r the set of ROCs to write
         * @param _qq_field the name of the qual field (this gets written into each row of the table)
         * @param _output_rocs true to write ROCs (otherwise, just write totals)
         * @param _roc_delta minimum QQ distance between ROC datapoints
         * @param _regions regions datastructure for adding region sizes
         */
        ROCOutput(RocMap const & r,
                  std::string _qq_field,
                  bool _output_rocs,
                  double _roc_delta,
                  variant::QuantifyRegions const & _regions,
                  std::vector<std::string> const & _roc_regions
        ) : rocs(r),
            qq_field(_qq_field),
            output_rocs(_output_rocs),
            roc_delta(_roc_delta),
            regions(_regions),
            roc_regions(_roc_regions)
        {}

        /** compute metrics and write output */
        void write(std::ostream & o);
    private:
        RocMap rocs;
        std::string qq_field;
        bool output_rocs;
        double roc_delta;
        variant::QuantifyRegions const & regions;
        std::vector<std::string> const & roc_regions;
   };
}

#endif //HAPLOTYPES_ROCOUTPUT_HH
