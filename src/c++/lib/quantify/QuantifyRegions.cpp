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
 * Track named regions
 *
 * \file QuantifyRegions.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/filesystem.hpp>

#include "QuantifyRegions.hh"

#include "helpers/IntervalBuffer.hh"
#include "helpers/BCFHelpers.hh"

#include <map>
#include <unordered_map>
#include <htslib/vcf.h>

#include "Fasta.hh"
#include "RefVar.hh"
#include "Error.hh"

namespace variant
{
    struct QuantifyRegions::QuantifyRegionsImpl
    {
        QuantifyRegionsImpl(std::string const &_ref) : ref(_ref.c_str())
        {}

        std::vector<std::string> names;
        std::vector<size_t> levels;
        std::unordered_map<std::string, size_t> label_map;
        std::unordered_map<std::string, std::unique_ptr<intervals::IntervalBuffer>> ib;
        std::unordered_map<std::string, std::unique_ptr<intervals::IntervalBuffer>>::iterator current_chr = ib.end();

        std::unordered_map<size_t, size_t> region_sizes;

        std::map< std::string, std::unordered_map<std::string, size_t> > extra_counts;

        int64_t current_pos = -1;
        FastaFile ref;
    };

    QuantifyRegions::QuantifyRegions(std::string const &ref) : _impl(new QuantifyRegionsImpl(ref))
    {
    }

    QuantifyRegions::~QuantifyRegions()
    {}

    /**
     * Returns true if regions were loaded.
     *
     * Note that this will also return true when an empty bed file was loaded.
     * This is intentional to distinguish the case where we don't have confident
     * regions (everything unknown is a FP) from the one where the confident
     * region file is empty (every FP is unknown).
     */
    bool QuantifyRegions::hasRegions(std::string const &rname) const
    {
        return _impl->label_map.find(rname) != _impl->label_map.cend();
    }

    /**
     * @return all region names
     */
    std::list<std::string> QuantifyRegions::getRegionNames() const
    {
        std::list<std::string> names {_impl->names.begin(), _impl->names.end()};
        return names;
    }

    void QuantifyRegions::load(std::vector<std::string> const &rnames, bool fixchr)
    {
        std::unordered_map<std::string, size_t> label_map;
        for (std::string const &f : rnames)
        {
            std::vector<std::string> v;
            stringutil::split(f, v, ":");
            bool fixed_label = false;

            std::string filename, label = "";

            // in case someone passes a ":"
            if (v.size() == 0)
            {
                error("Invalid region name: %s", f.c_str());
            }

            if (v.size() > 1)
            {
                label = v[0];
                filename = v[1];
                if (label[0] == '=')
                {
                    label = label.substr(1);
                    fixed_label = true;
                }
                if (label.find("CONF") == 0)
                {
                    // squash all CONF regions into one
                    label = "CONF";
                    fixed_label = true;
                }
            }
            else
            {
                filename = v[0];
                label = boost::filesystem::path(filename).stem().string();
            }

            htsFile *bedfile = NULL;

            if (stringutil::endsWith(filename, ".gz"))
            {
                bedfile = hts_open(filename.c_str(), "rz");
            }
            else
            {
                bedfile = hts_open(filename.c_str(), "r");
            }

            size_t icount = 0;
            kstring_t l;
            l.l = l.m = 0;
            l.s = NULL;
            size_t label_id;
            auto li_it = label_map.find(label);
            if (li_it == label_map.end())
            {
                label_id = _impl->names.size();
                _impl->names.push_back(label);
                _impl->levels.push_back(0);
                label_map[label] = label_id;
            }
            else
            {
                label_id = li_it->second;
            }
            while (hts_getline(bedfile, 2, &l) > 0)
            {
                std::string line(l.s);
                v.clear();
                stringutil::split(line, v, "\t");
                // we want >= 3 columns
                if (v.size() >= 3)
                {
                    if (fixchr)
                    {
                        if (v[0].size() > 0 && (
                            v[0].at(0) == '1' ||
                            v[0].at(0) == '2' ||
                            v[0].at(0) == '3' ||
                            v[0].at(0) == '4' ||
                            v[0].at(0) == '5' ||
                            v[0].at(0) == '6' ||
                            v[0].at(0) == '7' ||
                            v[0].at(0) == '8' ||
                            v[0].at(0) == '9' ||
                            v[0].at(0) == 'X' ||
                            v[0].at(0) == 'Y' ||
                            v[0].at(0) == 'M'))
                        {
                            v[0] = std::string("chr") + v[0];
                        }
                    }
                    auto chr_it = _impl->ib.find(v[0]);
                    if (chr_it == _impl->ib.end())
                    {
                        chr_it = _impl->ib.emplace(
                            v[0],
                            std::move(
                                std::unique_ptr<intervals::IntervalBuffer>(new intervals::IntervalBuffer()))).first;
                    }
                    // intervals are both zero-based
                    try
                    {

                        int64_t start = (int64_t) std::stoll(v[1]), stop = (int64_t) (std::stoll(v[2]) - 1);
                        if (start > stop)
                        {
                            std::cerr << "[W] ignoring invalid interval in " << filename << " : " << line << "\n";
                            continue;
                        }

                        /*
                         * create or get label id
                         */
                        auto getLabelId = [&label_map, this](std::string const & label_name, size_t level) -> size_t
                        {
                            auto li_it2 = label_map.find(label_name);
                            size_t this_label_id = 0;
                            if (li_it2 == label_map.end())
                            {
                                this_label_id = _impl->names.size();
                                _impl->names.push_back(label_name);
                                _impl->levels.push_back(level);
                                label_map[label_name] = this_label_id;
                            }
                            else
                            {
                                this_label_id = li_it2->second;
                            }
                            return this_label_id;
                        };

                        std::set<size_t> label_ids = {label_id};

                        if (!fixed_label && v.size() > 3)
                        {
                            size_t split = v[3].size();
                            for (size_t pos = v[3].size(); pos != 0; --pos)
                            {
                                if (v[3][pos] < '0' || v[3][pos] > '9')
                                {
                                    break;
                                }
                                else
                                {
                                    split = pos;
                                }
                            }

                            if(split < v[3].size() && v[3][split] == '_')
                            {
                                label_ids.insert(getLabelId(label + "_" + v[3].substr(0, split), 1));
                                label_ids.insert(getLabelId(label + "_" + v[3], 2));
                            }
                            else
                            {
                                label_ids.insert(getLabelId(label + "_" + v[3], 1));
                            }
                        }

                        for(const auto this_label_id : label_ids)
                        {
                            auto size_it = _impl->region_sizes.find(this_label_id);
                            if (size_it == _impl->region_sizes.end())
                            {
                                _impl->region_sizes[this_label_id] = (unsigned long) (stop - start + 1);
                            }
                            else
                            {
                                size_it->second += (unsigned long) (stop - start + 1);
                            }
                            chr_it->second->addInterval(start, stop, this_label_id);

                        }

                        ++icount;
                    }
                    catch (std::invalid_argument const &)
                    {
                        std::cerr << "[W] ignoring invalid interval in " << filename << " : " << line << "\n";
                    }
                    catch (std::out_of_range const &)
                    {
                        std::cerr << "[W] ignoring invalid interval in " << filename << " : " << line << "\n";
                    }
                }
                else if (line != "" && line != "\n")
                {
                    std::cerr << "[W] ignoring mis-formatted input line in " << filename << " : " << line << "\n";
                }
            }
            free(l.s);
            hts_close(bedfile);
            std::cerr << "Added region file '" << filename << "' as '" << label << "' (" << icount << " intervals)" <<
                      "\n";
        }



        _impl->label_map = label_map;
    }

    /** add Regions annotation to a record
     *
     * Records must be passed in sorted order.
     *
     */
    void QuantifyRegions::annotate(bcf_hdr_t *hdr, bcf1_t *record)
    {
        const std::string chr = bcfhelpers::getChrom(hdr, record);
        int64_t refstart = 0, refend = 0;
        bcfhelpers::getLocation(hdr, record, refstart, refend);

        const std::string ref_allele = record->d.allele[0];

        // pure insertion records must be fully contained in CONF to match
        bool is_pure_insertion = false;

        if (bcfhelpers::classifyAlleleString(ref_allele).first == bcfhelpers::AlleleType::NUC)
        {
            is_pure_insertion = record->n_allele > 1;
            int64_t updated_ref_start = std::numeric_limits<int64_t>::max();
            int64_t updated_ref_end = std::numeric_limits<int64_t>::min();
            bool nuc_alleles = false;

            for (int al = 1; al < record->n_allele; ++al)
            {
                RefVar al_rv;
                al_rv.start = record->pos;
                al_rv.end = (int64_t) (record->pos + ref_allele.size() - 1);
                auto ca = bcfhelpers::classifyAlleleString(record->d.allele[al]);
                if (ca.first == bcfhelpers::AlleleType::MISSING)
                {
                    ca.second = "";
                }
                else if (ca.first != bcfhelpers::AlleleType::NUC)
                {
                    break;
                }
                nuc_alleles = true;
                al_rv.alt = ca.second;

                variant::trimRight(_impl->ref, chr.c_str(), al_rv, false);
                variant::trimLeft(_impl->ref, chr.c_str(), al_rv, false);

                if (al_rv.end >= al_rv.start)
                {
                    updated_ref_start = std::min(updated_ref_start, al_rv.start);
                    updated_ref_end = std::max(updated_ref_end, al_rv.end);
                    is_pure_insertion = false;
                }
                else
                {
                    // this is an insertion *before* start, it will have end < start
                    // insertions are captured by the reference bases before and after
                    updated_ref_start = std::min(updated_ref_start, al_rv.start - 1);
                    updated_ref_end = std::max(updated_ref_end, al_rv.start);
                }
            }

            if (nuc_alleles)
            {
                refstart = updated_ref_start;
                refend = updated_ref_end;
            }
        }

        std::string tag_string = "";
        std::set<std::string> regions;

        const std::string regions_already_in_place = bcfhelpers::getInfoString(hdr, record, "Regions", "");
        std::vector<std::string> regions_split;
        stringutil::split(regions_already_in_place, regions_split, ",");
        for(const auto & r : regions_split)
        {
            // these we replace here
            if(r != "CONF" && r != "TS_boundary")
            {
                regions.insert(r);
            }
        }

        auto p_chr = _impl->current_chr;
        if (p_chr == _impl->ib.end() || p_chr->first != chr)
        {
            _impl->current_pos = -1;
            p_chr = _impl->ib.find(chr);
        }

        if (p_chr != _impl->ib.end())
        {
            if (refstart < _impl->current_pos)
            {
                error("Variants out of order at %s:%i", chr.c_str(), refstart);
            }
            for (size_t i = 0; i < _impl->names.size(); ++i)
            {
                if (p_chr->second->hasOverlap(refstart, refend, i))
                {
                    const bool fully_covered = p_chr->second->isCovered(refstart, refend, i);
                    if(!is_pure_insertion || fully_covered)
                    {
                        regions.insert(_impl->names[i]);
                        if (_impl->names[i] != "CONF" && !fully_covered)
                        {
                            regions.insert("TS_boundary");
                        }
                    }
                }
            }
            if (refstart > 1)
            {
                _impl->current_pos = refstart - 1;
                p_chr->second->advance(refstart - 1);
            }
        }
        // regions set is sorted, make sure Regions is sorted also
        for (auto const &r : regions)
        {
            if (!tag_string.empty())
            {
                tag_string += ",";
            }
            tag_string += r;
        }
        if (!tag_string.empty())
        {
            bcf_update_info_string(hdr, record, "Regions", tag_string.c_str());
        }
        else
        {
            bcf_update_info_string(hdr, record, "Regions", nullptr);
        }
        bcf_update_info_string(hdr, record, "RegionsExtent", (std::to_string(refstart + 1) + "-" +
                                                              std::to_string(refend + 1)).c_str());
    }

    /**
     * Get total region sizes in NT
     * @param region_name
     * @return  the region size
     */
    size_t QuantifyRegions::getRegionSize(std::string const &region_name) const
    {
        if(region_name == "*" || region_name == "TS_boundary")
        {
            return _impl->ref.contigNonNSize();
        }
        auto label_it = _impl->label_map.find(region_name);
        if(region_name == "TS_contained")
        {
            label_it = _impl->label_map.find("CONF");
        }

        if (label_it == _impl->label_map.cend())
        {
            return 0;
        }
        auto size_it = _impl->region_sizes.find(label_it->second);
        if (size_it == _impl->region_sizes.cend())
        {
            return 0;
        }
        return size_it->second;
    }

    /**
     * Get level for a region (0 for top-level, 1 for regions specified in a bed file)
     * @param region_name
     * @return level for the given region
     */
    size_t QuantifyRegions::getRegionLevel(std::string const & region_name) const
    {
        auto label = _impl->label_map.find(region_name);

        if(label == _impl->label_map.end())
        {
            return 0;
        }

        return _impl->levels[label->second];
    }

    /**
     * Mark one region to intersect with all others. This will
     * @param base base region (e.g. "CONF")
     */
    void QuantifyRegions::setIntersectRegion(std::string const & base)
    {
        auto base_label = _impl->label_map.find(base);

        if(base_label == _impl->label_map.end())
        {
            error("Unknown region label %s", base.c_str());
        }

        _impl->extra_counts[std::string("Subset.IS_") + base + ".Size"] = getRegionIntersectionSize(base);
        _impl->extra_counts[std::string("Subset.IS_") + base + ".Size"]["*"] = getRegionSize("CONF");
    }

    /**
     * Get all intersection / extra counts for a given region. By default, will return one
     * pair with <"Subset.Size", size of the region>
     * @param r the region to query counts for
     * @return a list of counts and names
     */
    std::list< std::pair<std::string, size_t> > QuantifyRegions::getRegionExtraCounts(std::string const & r) const
    {
        std::list< std::pair<std::string, size_t> > result;

        result.emplace_back("Subset.Size", getRegionSize(r));

        for(auto const & ec : _impl->extra_counts)
        {
            auto val = ec.second.find(r);
            if(val != ec.second.end())
            {
                result.emplace_back(ec.first, val->second);
            }
        }

        result.emplace_back("Subset.Level", getRegionLevel(r));

        return result;
    }

    /**
     * Get region intersection sizes in NT
     * @param region_name base region to intersect all others with (e.g. "CONF")
     * @return  mapping of sizes and overlaps by region name
     */
    std::unordered_map<std::string, size_t> QuantifyRegions::getRegionIntersectionSize(std::string const & region_name) const
    {
        std::unordered_map<std::string, size_t> result;

        auto this_label = _impl->label_map.find(region_name);
        if(this_label == _impl->label_map.end())
        {
            return result;
        }

        for(auto const & l : _impl->label_map)
        {
            if(l.first == region_name)
            {
                result[l.first] = _impl->region_sizes[l.second];
            }
            else
            {
                size_t total = 0;

                for(auto const & chr : _impl->ib)
                {
                    total += chr.second->intersectLanes(this_label->second, l.second);
                }
                result[l.first] = total;
            }
        }

        result["TS_boundary"] = result["CONF"];
        result["TS_contained"] = result["CONF"];

        return result;
    }
}


