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
 * Count variants
 *
 * \file BlockQuantify.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#include "helpers/StringUtil.hh"
#include "helpers/BCFHelpers.hh"
#include "Variant.hh"
#include "Fasta.hh"
#include "Error.hh"

#include <htslib/hts.h>

#include <array>
#include <list>
#include <htslib/vcf.h>
#include <thread>
#include <cmath>

#include "BlockQuantify.hh"
#include "BlockQuantifyImpl.hh"

#include "XCmpQuantify.hh"
#include "GA4GHQuantify.hh"

//#define DEBUG_BLOCKQUANTIFY

namespace variant {

    /** factory method */
    std::unique_ptr<BlockQuantify> makeQuantifier(
        bcf_hdr_t * hdr,
        FastaFile const & ref_fasta,
        std::string const & type,
        std::string const & params)
    {
        if(type == "xcmp")
        {
            return std::unique_ptr<BlockQuantify>((BlockQuantify*)new XCMPQuantify(hdr, ref_fasta, params));
        }
        else if(type == "ga4gh")
        {
            return std::unique_ptr<BlockQuantify>((BlockQuantify*)new GA4GHQuantify(hdr, ref_fasta, params));
        }
        else
        {
            error("Unknown quantification method: '%s'", type.c_str());
        }
        return nullptr;
    }

    /** list all quantification methods */
    std::list<std::string> listQuantificationMethods()
    {
        return {"xcmp", "ga4gh"};
    }

    BlockQuantify::BlockQuantifyImpl::~BlockQuantifyImpl()
    {
        for(auto x : variants)
        {
            bcf_destroy(x);
        }
    }

    BlockQuantify::BlockQuantify(bcf_hdr_t * hdr,
                                 FastaFile const & ref_fasta,
                                 std::string const & params) :
        _impl(std::unique_ptr<BlockQuantifyImpl>(new BlockQuantifyImpl{hdr,
                                                    ref_fasta,
                                                    nullptr,
                                                    BlockQuantifyImpl::count_map_t(),
                                                    BlockQuantifyImpl::variantlist_t(),
                                                    bcfhelpers::getSampleNames(hdr),
                                                    BlockQuantifyImpl::rocmap_t(),
                                                    BlockQuantifyImpl::filterset_t(),
                                                    params.find("count_unk") != std::string::npos,
                                                    params.find("output_vtc") != std::string::npos,
                                                    params.find("count_homref") != std::string::npos,
                                                    params.find("extended_counts") != std::string::npos
        }))
    {
    }

    BlockQuantify::~BlockQuantify() {}

    BlockQuantify::BlockQuantify(BlockQuantify && rhs) : _impl(std::move(rhs._impl)) { }

    BlockQuantify & BlockQuantify::operator=(BlockQuantify && rhs)
    {
        _impl = std::move(rhs._impl);
        return *this;
    }

    void BlockQuantify::rocFiltering(std::string const &filters)
    {
        std::vector<std::string> _fs;
        stringutil::split(filters, _fs, " ;,");
        _impl->filters_to_ignore.insert(_fs.cbegin(), _fs.cend());
    }

    // add ROC decision point
    void BlockQuantify::addROCValue(std::string const & roc_identifier,
                                    roc::DecisionType dt,
                                    double level,
                                    uint64_t n,
                                    bcf1_t * v)
    {
        // add observation to a roc
        auto observe = [this, level, dt, n, v](std::string const & name, bool f) {
            roc::DecisionType final_dt = dt;
            if(f)
            {
                if(dt == roc::DecisionType::TP)
                {
                    // filter-failed TPs become FNs
                    final_dt = roc::DecisionType::FN;
                }
                else if(dt == roc::DecisionType::TP2)
                {
                    // filter-failed TPs become FNs
                    final_dt = roc::DecisionType::FN2;
                }
                else if(dt != roc::DecisionType::FN)
                {
                    // filter-failed FPs / UNKs become Ns
                    final_dt = roc::DecisionType::N;
                }
            }

            uint64_t flags = 0;

            switch(final_dt)
            {
                case roc::DecisionType::FN:
                case roc::DecisionType::TP:
                {
                    const std::string bk_truth = bcfhelpers::getFormatString(_impl->hdr, v, "BK", 0, ".");
                    const std::string bi_truth = bcfhelpers::getFormatString(_impl->hdr, v, "BI", 0, ".");
                    const std::string blt_truth = bcfhelpers::getFormatString(_impl->hdr, v, "BLT", 0, ".");
                    flags = roc::makeObservationFlags(bk_truth, bi_truth, blt_truth);
                    break;
                }
                case roc::DecisionType::FP:
                case roc::DecisionType::UNK:
                case roc::DecisionType::TP2:
                case roc::DecisionType::FN2:
                {
                    const std::string bk_query = bcfhelpers::getFormatString(_impl->hdr, v, "BK", 1, ".");
                    const std::string bi_query = bcfhelpers::getFormatString(_impl->hdr, v, "BI", 1, ".");
                    const std::string blt_query = bcfhelpers::getFormatString(_impl->hdr, v, "BLT", 1, ".");
                    flags = roc::makeObservationFlags(bk_query, bi_query, blt_query);
                    break;
                }
                default: break;
            }

            auto it = _impl->rocs.find(name);
            if(it == _impl->rocs.end())
            {
                it = _impl->rocs.insert(std::make_pair(name, roc::Roc())).first;
            }

            // make sure FN and N always come first
            if( (   final_dt == roc::DecisionType::FN
                 || final_dt == roc::DecisionType::FN2
                 || final_dt == roc::DecisionType::N )
                && level == 0
                )
            {
                it->second.add(roc::Observation{std::numeric_limits<double>::min(), final_dt, n, flags});
            }
            else
            {
                it->second.add(roc::Observation{level, final_dt, n, flags});
            }
        };

        bcf_unpack(v, BCF_UN_FLT);
        bool fail = false;  // fails any of the non-blocked filters
        bool fail_any = false;  // fails any filter
        for(int j = 0; j < v->d.n_flt; ++j)
        {
            std::string filter = "PASS";
            int k = v->d.flt[j];
            if(k >= 0)
            {
                filter = bcf_hdr_int2id(_impl->hdr, BCF_DT_ID, v->d.flt[j]);
            }
            if(filter != "PASS" && filter != "")
            {
                if(!_impl->filters_to_ignore.count("*") && !_impl->filters_to_ignore.count(filter))
                {
                    fail = true;
                    observe("f:" + roc_identifier + ":" + filter, false);
                }
                else
                {
                    observe("f:" + roc_identifier + ":SEL_IGN_" + filter, false);
                }
                fail_any = true;
            }
        }

        observe("a:" + roc_identifier + ":PASS", fail_any);
        // selectively-filtered ROCs
        if(!_impl->filters_to_ignore.empty())
        {
            observe("a:" + roc_identifier + ":SEL", fail);
        }
        observe("a:" + roc_identifier + ":ALL", false);

        const std::string regions = bcfhelpers::getInfoString(_impl->hdr, v, "Regions", "");
        if(!regions.empty() && regions != "CONF")
        {
            std::vector<std::string> rs;
            stringutil::split(regions, rs, ",");
            for(auto const & r : rs)
            {
                if(r == "CONF")
                {
                    continue;
                }
                observe("s|" + r + ":" + roc_identifier + ":PASS", fail_any);
                if(!_impl->filters_to_ignore.empty())
                {
                    observe("s|" + r + ":" + roc_identifier + ":SEL", fail);
                }
                observe("s|" + r + ":" + roc_identifier + ":ALL", false);
            }
        }
    }

    void BlockQuantify::add(bcf1_t * v)
    {
        _impl->variants.push_back(v);
    }

    void BlockQuantify::count()
    {
        _impl->fasta_to_use.reset(new FastaFile(_impl->ref_fasta));
#ifdef DEBUG_BLOCKQUANTIFY
        int lastpos = 0;
        std::cerr << "starting block." << "\n";
#endif
        auto current_bs_start = _impl->variants.begin();
        std::string current_chr;
        int current_bs = -1;
        bool current_bs_valid = false;

        // function to compute the QQ values for truth variants in the current
        // benchmarking superlocus
        const auto update_bs_filters = [this, &current_bs_start](BlockQuantifyImpl::variantlist_t::iterator to)
        {
            std::set<int> bs_filters;
            for(auto cur = current_bs_start; cur != to; ++cur)
            {
                for(int nf = 0; nf < (*cur)->d.n_flt; ++nf)
                {
                    const int f = (*cur)->d.flt[nf];
                    if(f != bcf_hdr_id2int(_impl->hdr, BCF_DT_ID, "PASS"))
                    {
                        bs_filters.insert(f);
                    }
                }
            }

            if(bs_filters.empty())
            {
                return;
            }

            for(auto cur = current_bs_start; cur != to; ++cur)
            {
                const std::string bdt = bcfhelpers::getFormatString(_impl->hdr, *cur, "BD", 0);
                const std::string bvq = bcfhelpers::getFormatString(_impl->hdr, *cur, "BVT", 1);
                // filter TPs where the query call in NOCALL
                if(bdt == "TP" && bvq == "NOCALL")
                {
                    for(auto f : bs_filters)
                    {
                        bcf_add_filter(_impl->hdr, *cur, f);
                    }
                }
            }
        };

        // function to compute the QQ values for truth variants in the current
        // benchmarking superlocus
        const auto update_bs_qq = [this, &current_bs_start](BlockQuantifyImpl::variantlist_t::iterator to)
        {
            std::vector<float> tp_qqs;
            for(auto cur = current_bs_start; cur != to; ++cur)
            {
                const float qqq = bcfhelpers::getFormatFloat(_impl->hdr, *cur, "QQ", 1);
                if(std::isnan(qqq))
                {
                    continue;
                }
                const std::string bd = bcfhelpers::getFormatString(_impl->hdr, *cur, "BD", 1);
                // we want the scores of all TPs in this BS
                if(bd == "TP")
                {
                    tp_qqs.push_back(qqq);
                }
            }

            float t_qq = bcfhelpers::missing_float();
            if(!tp_qqs.empty())
            {
                t_qq = *(std::min_element(tp_qqs.begin(), tp_qqs.end()));
            }

            /** compute the median over all variants */
            int fsize = bcf_hdr_nsamples(_impl->hdr);
            float * fmt = (float*)calloc((size_t) fsize, sizeof(float));
            for(auto cur = current_bs_start; cur != to; ++cur)
            {
                const std::string bd = bcfhelpers::getFormatString(_impl->hdr, *cur, "BD", 0);
                bcf_get_format_float(_impl->hdr, *cur, "QQ", &fmt, &fsize);
                if(bd != "TP")
                {
                    fmt[0] = bcfhelpers::missing_float();
                }
                else
                {
                    const float qqq = bcfhelpers::getFormatFloat(_impl->hdr, *cur, "QQ", 1);
                    const std::string bd = bcfhelpers::getFormatString(_impl->hdr, *cur, "BD", 1);
                    if(bd == "TP" && !std::isnan(qqq))
                    {
                        fmt[0] = qqq;
                    }
                    else
                    {
                        fmt[0] = t_qq;
                    }

                }
                bcf_update_format_float(_impl->hdr, *cur, "QQ", fmt, fsize);
            }
            free(fmt);

#ifdef DEBUG_BLOCKQUANTIFY
            const int bs = bcfhelpers::getInfoInt(_impl->hdr, *current_bs_start, "BS", -1);
            std::string values;
            for(float x : tp_qqs)
            {
                values += std::to_string(x) + ",";
            }
            std::cerr << "BS: " << bs << " T_QQ = " << t_qq << " [" << values << "]" << "\n";
#endif
        };

        const auto update_bs_conf_boundary_flag = [this, &current_bs_start](BlockQuantifyImpl::variantlist_t::iterator to)
        {
            static const int has_conf = 1;
            static const int has_non_conf = 2;
            int conf_non_conf = 0;
            for(auto cur = current_bs_start; cur != to; ++cur)
            {
                const std::string regions = bcfhelpers::getInfoString(_impl->hdr, *cur, "Regions", "");

                if(regions.find("CONF") == std::string::npos)
                {
                    conf_non_conf |= has_non_conf;
                }
                else
                {
                    conf_non_conf |= has_conf;
                }
                if(regions.find("TS_boundary") != std::string::npos)
                {
                    conf_non_conf |= has_non_conf | has_conf;
                }
            }

            for(auto cur = current_bs_start; cur != to; ++cur)
            {
                const std::string regions = bcfhelpers::getInfoString(_impl->hdr, *cur, "Regions", "");

                if(conf_non_conf == (has_conf | has_non_conf))
                {
                    if(regions.find("TS_boundary") == std::string::npos)
                    {
                        bcf_update_info_string(_impl->hdr,
                                               *cur, "Regions",
                                               (regions.empty() ? "TS_boundary" :
                                                regions + ",TS_boundary").c_str());
                    }
                }
                else if(conf_non_conf == has_conf)
                {
                    if(regions.find("TS_contained") == std::string::npos)
                    {
                        // also flag fully confident superloci
                        bcf_update_info_string(_impl->hdr,
                                               *cur, "Regions",
                                               (regions.empty() ? "TS_contained" :
                                                regions + ",TS_contained").c_str());
                    }
                }
            }
        };


        for(auto v_it = _impl->variants.begin(); v_it != _impl->variants.end(); ++v_it)
        {
            // update fields, must output GA4GH-compliant fields
            countVariants(*v_it);

            // determine benchmarking superlocus
            const std::string vchr = bcfhelpers::getChrom(_impl->hdr, *v_it);
            const int vbs = bcfhelpers::getInfoInt(_impl->hdr, *v_it, "BS");
            if(!current_bs_valid)
            {
                current_bs = vbs;
                current_chr = vchr;
                current_bs_valid = true;
            }

#ifdef DEBUG_BLOCKQUANTIFY
            std::cerr << "current BS = " << current_bs << " vbs = " << vbs << "\n";
#endif

            if(   current_bs_start != v_it
               && (vbs != current_bs || vbs < 0 || vchr != current_chr))
            {
#ifdef DEBUG_BLOCKQUANTIFY
                std::cerr << "finishing BS = " << current_bs << " vbs = " << vbs << "\n";
#endif
                update_bs_qq(v_it);
                update_bs_filters(v_it);
                update_bs_conf_boundary_flag(v_it);
                current_bs = vbs;
                current_chr = vchr;
                current_bs_start = v_it;
            }
        }

        // do final superlocus (if any)
        update_bs_qq(_impl->variants.end());
        update_bs_filters(_impl->variants.end());
        update_bs_conf_boundary_flag(_impl->variants.end());

        for(auto & v : _impl->variants)
        {
#ifdef DEBUG_BLOCKQUANTIFY
            lastpos = v->pos;
#endif
            // use BD and BVT to make ROCs
            rocEvaluate(v);
        }
#ifdef DEBUG_BLOCKQUANTIFY
        std::cerr << "finished block " << lastpos << " - " << _impl->variants.size() << " records on thread " << std::this_thread::get_id() << "\n";
#endif
        _impl->fasta_to_use.reset(nullptr);
    }

    // GA4GH-VCF field-based ROC counting
    void BlockQuantify::rocEvaluate(bcf1_t * v)
    {
        if (_impl->samples.size() != 2)
        {
            // number of samples must be two, first one is truth, second is query
            return;
        }
        const std::string bd_truth = bcfhelpers::getFormatString(_impl->hdr, v, "BD", 0, ".");
        const std::string bd_query = bcfhelpers::getFormatString(_impl->hdr, v, "BD", 1, ".");
        const std::string vt_truth = bcfhelpers::getFormatString(_impl->hdr, v, "BVT", 0, ".");
        const std::string vt_query = bcfhelpers::getFormatString(_impl->hdr, v, "BVT", 1, ".");

        if(vt_truth != "NOCALL")
        {
            double qq = bcfhelpers::getFormatFloat(_impl->hdr, v, "QQ", 0);
            if(std::isnan(qq))
            {
                qq = 0;
            }

            if(bd_truth == "TP")
            {
                addROCValue(vt_truth, roc::DecisionType::TP, qq, 1, v);
            }
            else if(bd_truth == "FN")
            {
                addROCValue(vt_truth, roc::DecisionType::FN, qq, 1, v);
            }
        }

        if(vt_query != "NOCALL")
        {
            double qq = bcfhelpers::getFormatFloat(_impl->hdr, v, "QQ", 1);
            if(std::isnan(qq))
            {
                qq = 0;
            }
            if(bd_query == "FP")
            {
                addROCValue(vt_query, roc::DecisionType::FP, qq, 1, v);
            }
            else if(bd_query == "TP")
            {
                addROCValue(vt_query, roc::DecisionType::TP2, qq, 1, v);
            }
            else if(bd_query == "UNK")
            {
                addROCValue(vt_query, roc::DecisionType::UNK, qq, 1, v);
            }
        }
    }

    // result output
    std::map<std::string, VariantStatistics> const & BlockQuantify::getCounts() const
    {
        return _impl->count_map;
    }

    std::map<std::string, roc::Roc> const & BlockQuantify::getRocs() const
    {
        return _impl->rocs;
    }

    std::list<bcf1_t*> const & BlockQuantify::getVariants()
    {
        return _impl->variants;
    }

}
