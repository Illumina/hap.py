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


#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "helpers/StringUtil.hh"
#include "Variant.hh"
#include "Error.hh"

#include <htslib/hts.h>

#include <array>
#include <list>

#include "BlockQuantify.hh"

namespace variant {

    void QuantifyRegions::load(std::vector<std::string> const &rnames) {
        std::map<std::string, std::vector<Interval<std::string, int64_t> > > intervals;
        for (std::string const &f : rnames) {
            std::vector<std::string> v;
            stringutil::split(f, v, ":");

            std::string filename, label = "";

            // in case someone passes a ":"
            if (v.size() == 0) {
                error("Invalid region name: %s", f.c_str());
            }

            if (v.size() > 1) {
                label = v[0];
                filename = v[1];
            }
            else {
                filename = v[0];
                label = boost::filesystem::path(filename).stem().string();
            }

            htsFile *bedfile = NULL;

            if (stringutil::endsWith(filename, ".gz")) {
                bedfile = hts_open(filename.c_str(), "rz");
            }
            else {
                bedfile = hts_open(filename.c_str(), "r");
            }

            size_t icount = 0;
            kstring_t l;
            l.l = l.m = 0;
            l.s = NULL;
            while (hts_getline(bedfile, 2, &l) > 0) {
                std::string line(l.s);
                v.clear();
                stringutil::split(line, v, "\t");
                // we want >= 3 columns
                if (v.size() >= 3) {
                    auto chr_it = intervals.find(v[0]);
                    if (chr_it == intervals.end()) {
                        chr_it = intervals.emplace(
                            v[0],
                            std::vector<Interval<std::string, int64_t> >()).first;
                    }
                    // intervals are both zero-based
                    int64_t start = (int64_t) std::stoll(v[1]), stop = (int64_t) (std::stoll(v[2]) - 1);
                    if (start > stop) {
                        std::cerr << "[W] ignoring invalid interval in " << filename << " : " << line << "\n";
                    }
                    chr_it->second.push_back(Interval<std::string, int64_t>(start, stop, label));
                    ++icount;
                }
                else if (line != "" && line != "\n") {
                    std::cerr << "[W] ignoring mis-formatted input line in " << filename << " : " << line << "\n";
                }
            }
            free(l.s);
            hts_close(bedfile);
            std::cerr << "Added region file '" << filename << "' as '" << label << "' (" << icount << " intervals)" <<
            "\n";
        }

        for (auto &p : intervals) {
            regions[p.first] = IntervalTree<std::string, int64_t>(p.second);
        }
    }

    struct BlockQuantify::BlockQuantifyImpl {
        VariantReader & r;
        QuantifyRegions const & qregions;
        std::string ref_fasta;
        bool output_vtc;
        bool count_homref;

        typedef std::list< std::pair<std::string, std::string> > samplenames_t;
        typedef std::map<std::string, VariantStatistics> count_map_t;
        typedef std::list<Variants> variantlist_t;

        count_map_t count_map;
        variantlist_t variants;
        QuantifyRegions::regionmap_t::const_iterator current_chr;
        samplenames_t samples;
    };

    BlockQuantify::BlockQuantify(VariantReader & r,
                                 QuantifyRegions const & qregions,
                                 std::string const & ref_fasta,
                                 bool output_vtc,
                                 bool count_homref) :
        _impl(std::unique_ptr<BlockQuantifyImpl>(new BlockQuantifyImpl{r,
                                                    qregions,
                                                    ref_fasta,
                                                    output_vtc,
                                                    count_homref,
                                                    BlockQuantifyImpl::count_map_t(),
                                                    BlockQuantifyImpl::variantlist_t(),
                                                    qregions.regions.cend(),
                                                    BlockQuantifyImpl::samplenames_t()}))
    {
        r.getSampleList(_impl->samples);
    }

    BlockQuantify::~BlockQuantify() {}

    BlockQuantify::BlockQuantify(BlockQuantify && rhs) : _impl(std::move(rhs._impl)) { }

    BlockQuantify & BlockQuantify::operator=(BlockQuantify && rhs)
    {
        _impl = std::move(rhs._impl);
        return *this;
    }

    void BlockQuantify::add(Variants const &v)
    {
        _impl->variants.push_back(v);
    }

    void BlockQuantify::count()
    {
        for(auto & v : _impl->variants)
        {
            count_variants(v);
        }
    }

    // result output
    std::map<std::string, VariantStatistics> const & BlockQuantify::getCounts() const
    {
        return _impl->count_map;
    }

    std::list<Variants> const & BlockQuantify::getVariants() const
    {
        return _impl->variants;
    }

    /** save lookups when reading from the same chromosome */
    void BlockQuantify::start_chr(std::string const & chr) {
        if(   _impl->current_chr == _impl->qregions.regions.cend()
           || chr != _impl->current_chr->first)
        {
            _impl->current_chr = _impl->qregions.regions.find(chr);
        }
    }

    void BlockQuantify::count_variants(Variants & v) {
        // resolve tags
        std::set<std::string> to_count;
        std::string tag_string = "";

        start_chr(v.chr);
        if (_impl->current_chr != _impl->qregions.regions.cend()) {
            std::vector<Interval<std::string, int64_t> > overlapping;
            _impl->current_chr->second.findOverlapping(v.pos, v.pos + v.len - 1, overlapping);
            for (auto const &iv : overlapping) {
                to_count.insert(iv.value);
            }
            for (std::string const &r : to_count) {
                if (tag_string != "") {
                    tag_string += ",";
                }
                tag_string += r;
            }
        }

        std::string type = ".";
        std::string kind = ".";
        std::string gtt1 = ".";
        std::string gtt2 = ".";
        bool hapmatch = false;

        std::vector<std::string> infs;
        stringutil::split(v.info, infs, ";");
        v.info = "";
        for (std::string &i : infs) {
            if (i == "IMPORT_FAIL") {
                // preserve import fail annotation
                v.info += i + ";";
                continue;
            }
            if (i.size() < 2) {
                continue;
            }
            std::vector<std::string> ifo;
            stringutil::split(i, ifo, "=");
            if (ifo[0] == "HapMatch") {
                hapmatch = true;
            } else if (ifo.size() < 2) {
                v.info += i + ";";
            } else if (ifo[0] == "type") {
                type = ifo[1];
            } else if (ifo[0] == "kind") {
                kind = ifo[1];
            } else if (ifo[0] == "gtt1") {
                gtt1 = ifo[1];
            } else if (ifo[0] == "gtt2") {
                gtt2 = ifo[1];
            } else if (ifo[0] != "ctype") {
                v.info += i + ";";
            }
        }

        if (hapmatch && type != "TP") {
            kind = "hapmatch__" + type + "__" + kind;
            type = "TP";
        }

        if (!_impl->qregions.regions.empty() && tag_string.empty() && type == "FP" && kind == "missing") {
            type = "UNK";
        }

        std::string matchtypeinfo;
        matchtypeinfo += std::string("type=") + type;
        matchtypeinfo += std::string(";kind=") + kind;
        if (!tag_string.empty()) {
            matchtypeinfo += std::string(";Regions=") + tag_string;
        }

        std::set<int> vtypes;
        uint64_t hl_lt_truth = (uint64_t) -1;
        uint64_t hl_lt_query = (uint64_t) -1;
        uint64_t hl_vt_truth = (uint64_t) -1;
        uint64_t hl_vt_query = (uint64_t) -1;

        std::string lt_truth = "unk";
        std::string lt_query = "unk";

        // count all and specific type instance in
        // all samples
        const std::array<std::string, 2> names {"all", type + ":" + kind + ":" + tag_string};
        for (auto const & name : names)
        {
            int i = 0;
            for (auto const &s : _impl->samples) {
                std::string key;
                if (s.second == "") {
                    key = name + ":" + s.first + ":" + s.second + ":" + std::to_string(i);
                }
                else {
                    key = name + ":" + s.second;
                }

                auto it = _impl->count_map.find(key);
                if (it == _impl->count_map.end()) {
                    it = _impl->count_map.emplace(key, VariantStatistics(_impl->ref_fasta.c_str(),
                                                                         _impl->count_homref)).first;
                }
                int *types;
                int ntypes = 0;
                it->second.add(v, i, &types, &ntypes);
                uint64_t vts_seen = 0;
                uint64_t lt = (uint64_t) -1;
                for (int j = 0; j < ntypes; ++j) {
                    vts_seen |= types[j] & 0xf;
                    vtypes.insert(types[j]);
                    if ((types[j] & 0xf0) > 0x10) {
                        lt = (uint64_t) (types[j] >> 4);
                    }
                }
                if (s.second == "TRUTH") {
                    hl_vt_truth = vts_seen == VT_NOCALL ? 2 : ((vts_seen == variant::VT_SNP
                                                                || vts_seen == (variant::VT_SNP | variant::VT_REF)) ? 0
                                                                                                                    : 1);
                    hl_lt_truth = lt;
                } else if (s.second == "QUERY") {
                    hl_vt_query = vts_seen == VT_NOCALL ? 2 : ((vts_seen == variant::VT_SNP
                                                                || vts_seen == (variant::VT_SNP | variant::VT_REF)) ? 0
                                                                                                                    : 1);
                    hl_lt_query = lt;
                }
                ++i;
            }
        }

        static const char *nvs[] = {"UNK", "SNP", "INDEL", "NOCALL"};
        std::string type_info;
        if (hl_vt_truth < 4) {
            if(type_info != "") { type_info += ";"; }
            type_info += std::string("T_VT=") + nvs[hl_vt_truth + 1];
        }
        if (hl_vt_query < 4) {
            if (type_info != "") { type_info += ";"; }
            type_info += std::string("Q_VT=") + nvs[hl_vt_query + 1];
        }
        if (hl_lt_truth < 0x10) {
            if (type_info != "") { type_info += ";"; }
            type_info += std::string(";T_LT=") + CT_NAMES[hl_lt_truth];
        }
        if (hl_lt_query < 0x10) {
            if (type_info != "") { type_info += ";"; }
            type_info += std::string(";Q_LT=") + CT_NAMES[hl_lt_query];
        }
        if (_impl->output_vtc && !vtypes.empty()) {
            std::string s = "VTC=";
            int j = 0;
            for (int t : vtypes) {
                if (j++ > 0) { s += ","; }
                s += VariantStatistics::type2string(t);
            }
            if (type_info != "") { type_info += ";"; }
            type_info += s;
        }
        if(!type_info.empty())
        {
            if(!v.info.empty())
            {
                v.info += ";";
            }
            v.info += type_info;
        }
        if(!matchtypeinfo.empty())
        {
            if(!v.info.empty())
            {
                v.info += ";";
            }
            v.info += matchtypeinfo;
        }
    }
}
