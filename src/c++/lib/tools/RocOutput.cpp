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
 * Advanced ROC computations (Ti/Tv, AuC, ...)
 *
 * \file RocExtended.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <list>
#include <array>
#include <cmath>
#include <set>
#include <map>
#include <limits>
#include <array>

#include "helpers/RocOutput.hh"
#include "helpers/Table.hh"
#include "helpers/StringUtil.hh"
#include "Error.hh"


namespace roc
{
    struct KEYS
    {
        enum _KEYS {
            Type = 0,
            Subtype,
            Genotype,
            Filter,
            Subset,
            QQ_Field,
            QQ,
            SIZE
        };
    };

    static const std::array<std::string, KEYS::SIZE> SKEYS {
        "Type",
        "Subtype",
        "Genotype",
        "Filter",
        "Subset",
        "QQ.Field",
        "QQ"
    };

    struct METRICS
    {
        enum  _METRICS
        {
            TP = 0,
            FN,
            FP,
            TP2,
            UNK,
            TRUTH_TOTAL,
            QUERY_TOTAL,
            FP_gt,
            FP_al,
            Recall,
            Precision,
            F1_Score,
            Frac_NA,
            SIZE
        };
    };

    static const std::array<std::string, METRICS::SIZE> SMETRICS {
        "TRUTH.TP",
        "TRUTH.FN",
        "QUERY.FP",
        "QUERY.TP",
        "QUERY.UNK",
        "TRUTH.TOTAL",
        "QUERY.TOTAL",
        "FP.gt",
        "FP.al",
        "METRIC.Recall",
        "METRIC.Precision",
        "METRIC.F1_Score",
        "METRIC.Frac_NA",
    };

    static const std::list<std::pair<std::string, uint64_t> > SNP_SUBTYPES
    {
        {"*", 0},
        {"ti",  roc::OBS_FLAG_TI},
        {"tv",  roc::OBS_FLAG_TV},
    };

    static const std::list<std::pair<std::string, uint64_t> > INDEL_SUBTYPES
    {
        {"*", 0},
        {"I1_5",  roc::OBS_FLAG_I1_5},
        {"I6_15",  roc::OBS_FLAG_I6_15},
        {"I16_PLUS",  roc::OBS_FLAG_I16_PLUS},
        {"D1_5",  roc::OBS_FLAG_D1_5},
        {"D6_15",  roc::OBS_FLAG_D6_15},
        {"D16_PLUS",  roc::OBS_FLAG_D16_PLUS},
        {"C1_5",  roc::OBS_FLAG_C1_5},
        {"C6_15",  roc::OBS_FLAG_C6_15},
        {"C16_PLUS",  roc::OBS_FLAG_C16_PLUS},
    };

    template<class _e>
    static const std::string _S(_e e);

    template<>
    const std::string _S(KEYS::_KEYS e)
    {
        return SKEYS[e];
    }

    template<>
    const std::string _S(METRICS::_METRICS e)
    {
        return SMETRICS[e];
    }

    /** count Ti/Tv and het/hom ratios for these counts */
    static const std::list<METRICS::_METRICS> SUBTYPED_COUNTS = {
        METRICS::TP,
        METRICS::FN,
        METRICS::FP,
        METRICS::TP2,
        METRICS::UNK,
        METRICS::TRUTH_TOTAL,
        METRICS::QUERY_TOTAL,
    };

    /** helper to write values from a single level */
    void _addLevel(const std::string & type,
                   const std::string & subtype,
                   const std::string & genotype,
                   const std::string & filter,
                   const std::string & subset,
                   const std::string & qq_field,
                   Level const & l,
                   table::Table & table,
                   bool counts_only,
                   variant::QuantifyRegions const & regions)
    {
        // don't write rows for ti / tv / genotypes, but rather show the counts inline
        if(subtype != "ti" && subtype != "tv" && genotype == "*")
        {
            std::string prefix;
            if(std::isnan(l.level))
            {
                prefix =  type + "\t" + subtype + "\t*\t" + filter + "\t" + subset + "\t*";
                table.set(prefix, _S(KEYS::QQ), "*");
            }
            else
            {
                prefix = type + "\t" + subtype + "\t*\t"
                         + filter + "\t" + subset + "\t" + std::to_string(l.level);
                table.set(prefix, _S(KEYS::QQ), l.level);
            }
            table.set(prefix, _S(KEYS::Type), type);
            table.set(prefix, _S(KEYS::Subtype), subtype);
            table.set(prefix, _S(KEYS::Genotype), genotype);
            table.set(prefix, _S(KEYS::Subset), subset);
            table.set(prefix, _S(KEYS::Filter), filter);
            table.set(prefix, _S(KEYS::QQ_Field), qq_field);

            table.set(prefix, _S(METRICS::TP), l.tp());
            table.set(prefix, _S(METRICS::TP2), l.tp2());
            table.set(prefix, _S(METRICS::FP), l.fp());
            table.set(prefix, _S(METRICS::UNK), l.unk());
            table.set(prefix, _S(METRICS::FP_al), l.fp_al());
            table.set(prefix, _S(METRICS::FP_gt), l.fp_gt());

            const auto extra_counts = regions.getRegionExtraCounts(subset);

            for(auto const & ec : extra_counts)
            {
                table.set(prefix, ec.first, ec.second);
            }

            if(!counts_only)
            {
                table.set(prefix, _S(METRICS::FN), l.fn());
                table.set(prefix, _S(METRICS::Recall), l.recall());
                table.set(prefix, _S(METRICS::Precision), l.precision());
                table.set(prefix, _S(METRICS::F1_Score), l.fScore());
                table.set(prefix, _S(METRICS::Frac_NA), l.na());
                table.set(prefix, _S(METRICS::TRUTH_TOTAL), l.totalTruth());
                table.set(prefix, _S(METRICS::QUERY_TOTAL), l.totalQuery());
            }
        } else {
            // save ti and tv counts for computing ratios later
            if(genotype == "*" && (subtype == "ti" || subtype == "tv"))
            {
                std::string subtype_agg_prefix;
                if(std::isnan(l.level))
                {
                    subtype_agg_prefix =  type + "\t*\t*\t" + filter + "\t" + subset + "\t*";
                }
                else
                {
                    subtype_agg_prefix = type + "\t*\t*\t"
                             + filter + "\t" + subset + "\t" + std::to_string(l.level);
                }
                table.set(subtype_agg_prefix, _S(METRICS::TP) + "." + subtype, l.tp());
                table.set(subtype_agg_prefix, _S(METRICS::TP2) + "." + subtype, l.tp2());
                table.set(subtype_agg_prefix, _S(METRICS::FP) + "." + subtype, l.fp());
                table.set(subtype_agg_prefix, _S(METRICS::UNK) + "." + subtype, l.unk());
                if(!counts_only)
                {
                    table.set(subtype_agg_prefix, _S(METRICS::FN) + "." + subtype, l.fn());
                    table.set(subtype_agg_prefix, _S(METRICS::TRUTH_TOTAL) + "." + subtype, l.totalTruth());
                    table.set(subtype_agg_prefix, _S(METRICS::QUERY_TOTAL) + "." + subtype, l.totalQuery());
                }
            }

            // save het / hom counts for ratios
            if(genotype != "*" && subtype != "ti" && subtype != "tv")
            {
                std::string genotype_agg_prefix;
                if(std::isnan(l.level))
                {
                    genotype_agg_prefix =  type + "\t" + subtype + "\t*\t" + filter + "\t" + subset + "\t*";
                }
                else
                {
                    genotype_agg_prefix = type + "\t" + subtype + "\t*\t" + "\t"
                             + filter + "\t" + subset + "\t" + std::to_string(l.level);
                }
                table.set(genotype_agg_prefix, _S(METRICS::TP) + "." + genotype, l.tp());
                table.set(genotype_agg_prefix, _S(METRICS::TP2) + "." + genotype, l.tp2());
                table.set(genotype_agg_prefix, _S(METRICS::FP) + "." + genotype, l.fp());
                table.set(genotype_agg_prefix, _S(METRICS::UNK) + "." + genotype, l.unk());
                if(!counts_only)
                {
                    table.set(genotype_agg_prefix, _S(METRICS::FN) + "." + genotype, l.fn());
                    table.set(genotype_agg_prefix, _S(METRICS::TRUTH_TOTAL) + "." + genotype, l.totalTruth());
                    table.set(genotype_agg_prefix, _S(METRICS::QUERY_TOTAL) + "." + genotype, l.totalQuery());
                }
            }
        }
    }

    void ROCOutput::write(std::ostream &out_roc)
    {
        /** populate output table */
        table::Table output_values;

        std::set<std::string> filters;
        const std::list<std::pair<std::string, uint64_t> > gts = {
            {"het",    roc::OBS_FLAG_HET},
            {"hetalt", roc::OBS_FLAG_HETALT},
            {"homalt", roc::OBS_FLAG_HOMALT},
            {"*",      0},
        };

        const auto region_names = regions.getRegionNames();

        // create dummy entries for regions we don't have data for
        for(auto & m : rocs)
        {
            std::vector<std::string> subs;
            stringutil::split(m.first, subs, ":");
            if (subs.empty())
            {
                // only happens when an empty string is passed
                continue;
            }
            std::string type = "unknown";
            std::string subset = "unknown";
            std::string filter = "unknown";

            if (subs.size() == 3 && subs[0] == "a")
            {
                type = subs[1];
                filter = subs[2];
                for(auto const & region : region_names)
                {
                    // don't include CONF regions
                    if(region.find("CONF") == 0)
                    {
                        continue;
                    }
                    const std::string rname = std::string("s|") + region + ":" + type + ":" + filter;
                    auto it = rocs.find(rname);
                    if(it == rocs.end())
                    {
                        rocs[rname] = Roc();
                    }
                }
            }
        }

        for(auto & m : rocs)
        {
            // filter ROCs just give numbers of TPs / FPs / UNKs filtered by
            // each filter at each level

            std::vector<std::string> subs;
            stringutil::split(m.first, subs, ":");
            if (subs.empty())
            {
                // only happens when an empty string is passed
                continue;
            }
            std::string type = "unknown";
            std::string subset = "unknown";
            std::string filter = "unknown";

            bool counts_only = false;
            if (subs.size() != 3)
            {
                type = m.first;
            }
            else
            {
                type = subs[1];
                filter = subs[2];
                if (subs[0] == "f")
                {
                    counts_only = true;
                    subset = "*";
                }
                else if (subs[0].substr(0, 2) == "s|")
                {
                    subset = subs[0].substr(2);
                }
                else if (subs[0] == "a")
                {
                    subset = "*";
                }
                else
                {
                    // unknown count, should not happen
                    continue;
                }
            }

            counts_only = counts_only &&
                          (std::find(roc_regions.begin(),
                                     roc_regions.end(),
                                     subset) != roc_regions.end());

            std::list<std::pair<std::pair<std::string, std::string>, uint64_t> > subtypes;
            std::list<std::pair<std::string, uint64_t> > const * sts = nullptr;

            if (type == "SNP")
            {
                sts = &SNP_SUBTYPES;
            }
            else if (type == "INDEL")
            {
                sts = &INDEL_SUBTYPES;
            }
            else
            {
                continue;
            }

            for (auto const &st : *sts)
            {
                for (auto const &gt : gts)
                {
                    subtypes.emplace_back(std::make_pair(st.first, gt.first), st.second | gt.second);
                }
            }

            for (auto const &st : subtypes)
            {
                const roc::Level l = m.second.getTotals(st.second);
                _addLevel(type, st.first.first, st.first.second, filter, subset, qq_field, l, output_values, counts_only, regions);

                // don't write ROCs for filters
                // we could make this optional, but not sure it is really necessary
                if ((!counts_only) && output_rocs)
                {
                    std::vector<roc::Level> levels;
                    m.second.getLevels(levels, roc_delta, st.second);
                    for (auto const &l2 : levels)
                    {
                        _addLevel(type, st.first.first, st.first.second, filter, subset, qq_field, l2, output_values, counts_only, regions);
                    }
                }
            }
        }
        output_values.dropRowsWithMissing(_S(KEYS::Type));
        out_roc << output_values;
    }
}
