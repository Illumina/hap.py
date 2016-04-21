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
#include <cmath>
#include <set>
#include <map>
#include <limits>

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
            TiTv_ratio,
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
        "METRIC.TiTv_ratio",
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

    /** helper to write values from a single level */
    void _addLevel(const std::string & type,
                   const std::string & subtype,
                   const std::string & genotype,
                   const std::string & filter,
                   const std::string & subset,
                   const std::string & qq_field,
                   Level const & l,
                   table::Table & table,
                   bool counts_only)
    {
        std::string prefix;
        if(std::isnan(l.level))
        {
            prefix =  type + "\t" + subtype + "\t" + genotype + "\t" + filter + "\t" + subset + "\t*";
            table.set(prefix, _S(KEYS::QQ), "*");
        }
        else
        {
            prefix = type + "\t" + subtype + "\t" + genotype + "\t"
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
    }

    void ROCOutput::write(std::ostream &out_roc) const
    {
        /** populate output table */
        table::Table output_values;

        std::set<std::string> types;
        std::set<std::string> filters;
        std::set<std::string> subsets {"*"};
        std::set<std::string> qqs {"*"};
        const std::list<std::pair<std::string, uint64_t> > gts = {
            {"het",    roc::OBS_FLAG_HET},
            {"hetalt", roc::OBS_FLAG_HETALT},
            {"homalt", roc::OBS_FLAG_HOMALT},
            {"*",      0},
        };
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
            types.insert(type);
            subsets.insert(subset);
            filters.insert(filter);

            std::list<std::pair<std::pair<std::string, std::string>, uint64_t> > subtypes;
            std::list<std::pair<std::string, uint64_t> > const * sts = nullptr;

            // all types, all genotypes

            // SNP: TI / TV
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
                _addLevel(type, st.first.first, st.first.second, filter, subset, qq_field, l, output_values, counts_only);

                // don't write ROCs for filters
                // we could make this optional, but not sure it is really necessary
                if ((!counts_only) && output_rocs)
                {
                    std::vector<roc::Level> levels;
                    m.second.getLevels(levels, roc_delta, st.second);
                    for (auto const &l2 : levels)
                    {
                        _addLevel(type, st.first.first, st.first.second, filter, subset, qq_field, l2, output_values, counts_only);
                        qqs.insert(std::to_string(l2.level));
                    }
                }
            }
        }

        /** compute Ti/Tv ratios */
        const std::list<METRICS::_METRICS> counts = {
            METRICS::TP,
            METRICS::FN,
            METRICS::FP,
            METRICS::TP2,
            METRICS::UNK,
            METRICS::TRUTH_TOTAL,
            METRICS::QUERY_TOTAL,
        };

        for(auto const & subset : subsets)
        {
            for(auto const & filter : filters)
            {
                for(auto const & gt : gts)
                {
                    for(auto const & qq : qqs)
                    {
                        const std::string rowid = std::string("SNP\t*\t") + gt.first + "\t"
                                                  + filter + "\t" + subset + "\t" + qq;
                        const std::string rowid_ti = std::string("SNP\tti\t") + gt.first + "\t"
                                                     + filter + "\t" + subset + "\t" + qq;
                        const std::string rowid_tv = std::string("SNP\ttv\t") + gt.first + "\t"
                                                     + filter + "\t" + subset + "\t" + qq;
                        if(output_values.hasRow(rowid) &&
                           output_values.hasRow(rowid_ti) &&
                           output_values.hasRow(rowid_tv))
                        {
                            for(auto c : counts)
                            {
                                double ti_count = output_values.getDouble(rowid_ti, _S(c), 0);
                                double tv_count = output_values.getDouble(rowid_tv, _S(c), 0);
                                if(fabs(tv_count) < std::numeric_limits<double>::epsilon()
                                    || std::isnan(ti_count) || std::isnan(tv_count))
                                {
                                    output_values.set(rowid, _S(c) + ".TiTv_ratio", ".");
                                }
                                else
                                {
                                    output_values.set(rowid, _S(c) + ".TiTv_ratio", ti_count / tv_count);
                                }
                            }
                        }
                    }
                }
            }
        }

        for(auto const & type : types)
        {
            for(auto const & subset : subsets)
            {
                for(auto const & filter : filters)
                {
                    std::list<std::pair<std::string, uint64_t> >const * sts = NULL;
                    if(type == "SNP")
                    {
                        sts = &SNP_SUBTYPES;
                    }
                    else if(type == "INDEL")
                    {
                        sts = &INDEL_SUBTYPES;
                    }
                    else
                    {
                        continue;
                    }
                    for(auto const & st : *sts)
                    {
                        for(auto const & qq : qqs)
                        {
                            const std::string rowid = type + "\t" + st.first + "\t*" + "\t"
                                                      + filter + "\t" + subset + "\t" + qq;
                            const std::string rowid_het = type + "\t" + st.first + "\thet" + "\t"
                                                          + filter + "\t" + subset + "\t" + qq;
                            const std::string rowid_homalt = type + "\t" + st.first + "\thomalt" + "\t"
                                                          + filter + "\t" + subset + "\t" + qq;
                            if(output_values.hasRow(rowid) &&
                               output_values.hasRow(rowid_het) &&
                               output_values.hasRow(rowid_homalt))
                            {
                                for(auto c : counts)
                                {
                                    double het_count = output_values.getDouble(rowid_het, _S(c), 0);
                                    double homalt_count = output_values.getDouble(rowid_homalt, _S(c), 0);
                                    if(fabs(homalt_count) < std::numeric_limits<double>::epsilon()
                                       || std::isnan(het_count) || std::isnan(homalt_count))
                                    {
                                        output_values.set(rowid, _S(c) + ".het_hom_ratio", ".");
                                    }
                                    else
                                    {
                                        output_values.set(rowid, _S(c) + ".het_hom_ratio", het_count / homalt_count);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        out_roc << output_values;
    }
}
