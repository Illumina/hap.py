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
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "QuantifyRegions.hh"

#include "helpers/IntervalBuffer.hh"
#include "helpers/BCFHelpers.hh"

#include <map>
#include <unordered_map>

#include "Error.hh"

namespace variant
{
    struct QuantifyRegions::QuantifyRegionsImpl
    {
        std::vector<std::string> names;
        std::unordered_map<std::string, std::unique_ptr<intervals::IntervalBuffer>> ib;
        std::unordered_map<std::string, std::unique_ptr<intervals::IntervalBuffer>>::iterator current_chr = ib.end();
        int64_t current_pos = -1;
    };

    QuantifyRegions::QuantifyRegions() : _impl(new QuantifyRegionsImpl())
    { }

    QuantifyRegions::~QuantifyRegions()
    { }

    /**
     * Returns true if regions were loaded.
     *
     * Note that this will also return true when an empty bed file was loaded.
     * This is intentional to distinguish the case where we don't have confident
     * regions (everything unknown is a FP) from the one where the confident
     * region file is empty (every FP is unknown).
     */
    bool QuantifyRegions::hasRegions() const
    {
        return !_impl->names.empty();
    }

    void QuantifyRegions::load(std::vector<std::string> const &rnames)
    {
        std::map<std::string, size_t> label_map;
        for (std::string const &f : rnames)
        {
            std::vector<std::string> v;
            stringutil::split(f, v, ":");

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
                    auto chr_it = _impl->ib.find(v[0]);
                    if (chr_it == _impl->ib.end())
                    {
                        chr_it = _impl->ib.emplace(
                            v[0],
                            std::move(
                                std::unique_ptr<intervals::IntervalBuffer>(new intervals::IntervalBuffer()))).first;
                    }
                    // intervals are both zero-based
                    int64_t start = (int64_t) std::stoll(v[1]), stop = (int64_t) (std::stoll(v[2]) - 1);
                    if (start > stop)
                    {
                        std::cerr << "[W] ignoring invalid interval in " << filename << " : " << line << "\n";
                    }
                    chr_it->second->addInterval(start, stop, label_id);
                    ++icount;
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
    }

    /** add Regions annotation to a record
     *
     * Records must be passed in sorted order.
     *
     */
    void QuantifyRegions::annotate(bcf_hdr_t * hdr, bcf1_t *record)
    {
        std::string chr = bcfhelpers::getChrom(hdr, record);
        int64_t refstart = 0, refend = 0;
        bcfhelpers::getLocation(hdr, record, refstart, refend);

        std::string tag_string = "";

        auto p_chr = _impl->current_chr;
        if(p_chr == _impl->ib.end() || p_chr->first != chr)
        {
            _impl->current_pos = -1;
            p_chr = _impl->ib.find(chr);
        }

        if(p_chr != _impl->ib.end())
        {
            if(refstart < _impl->current_pos)
            {
                error("Variants out of order at %s:%i", chr.c_str(), refstart);
            }
            for(size_t i = 0; i < _impl->names.size(); ++i)
            {
                if(p_chr->second->hasOverlap(refstart, refend, i))
                {
                    if(!tag_string.empty())
                    {
                        tag_string += ",";
                    }
                    tag_string += _impl->names[i];
                }
            }
            if(refstart > 1)
            {
                _impl->current_pos = refstart - 1;
                p_chr->second->advance(refstart-1);
            }
        }
        if(!tag_string.empty())
        {
            bcf_update_info_string(hdr, record, "Regions", tag_string.c_str());
        }
        else
        {
            bcf_update_info_string(hdr, record, "Regions", nullptr);
        }
    }
}


