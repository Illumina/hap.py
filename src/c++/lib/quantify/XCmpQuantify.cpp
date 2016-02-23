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
 * \file XCmpQuantify.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "XCmpQuantify.hh"
#include "BlockQuantifyImpl.hh"

namespace variant
{
    XCMPQuantify::XCMPQuantify(bcf_hdr_t * hdr,
                               FastaFile const & ref_fasta,
                               const std::string & params) :
        BlockQuantify(hdr, ref_fasta, params)
    {
    }

    void XCMPQuantify::updateHeader(bcf_hdr_t * hdr)
    {
        // we need these, otherwise all the update info calls will fail
        bcf_hdr_append(hdr, "##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">");
        bcf_hdr_append(hdr, "##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">");
        bcf_hdr_append(hdr, "##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
        bcf_hdr_append(hdr, "##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
        bcf_hdr_append(hdr, "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">");
        bcf_hdr_append(hdr, "##INFO=<ID=T_VT,Number=1,Type=String,Description=\"High-level variant type in truth (SNP|INDEL).\">");
        bcf_hdr_append(hdr, "##INFO=<ID=Q_VT,Number=1,Type=String,Description=\"High-level variant type in query (SNP|INDEL).\">");
        bcf_hdr_append(hdr, "##INFO=<ID=T_LT,Number=1,Type=String,Description=\"High-level location type in truth (het|hom|hetalt).\">");
        bcf_hdr_append(hdr, "##INFO=<ID=Q_LT,Number=1,Type=String,Description=\"High-level location type in query (het|hom|hetalt).\">");
        if(_impl->output_vtc)
        {
            bcf_hdr_append(hdr, "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">");
        }
        bcf_hdr_sync(hdr);
    }

    void XCMPQuantify::countVariants(bcf1_t * v)
    {
        bcf_unpack(v, BCF_UN_ALL);
        std::string tag_string = bcfhelpers::getInfoString(_impl->hdr, v, "Regions", "");
        std::string type = bcfhelpers::getInfoString(_impl->hdr, v, "type");
        std::string kind = bcfhelpers::getInfoString(_impl->hdr, v, "kind");
        bool hapmatch = bcfhelpers::getInfoFlag(_impl->hdr, v, "HapMatch");

        if (hapmatch && type != "TP") {
            kind = "hapmatch__" + type + "__" + kind;
            type = "TP";
        }

        if (_impl->count_unk && tag_string.empty() && type == "FP" && kind == "missing") {
            type = "UNK";
        }

        bcf_update_info_string(_impl->hdr, v, "type", type.c_str());
        bcf_update_info_string(_impl->hdr, v, "kind", kind.c_str());

        std::set<int> vtypes;
        uint64_t hl_lt_truth = (uint64_t) -1;
        uint64_t hl_lt_query = (uint64_t) -1;
        uint64_t hl_vt_truth = (uint64_t) -1;
        uint64_t hl_vt_query = (uint64_t) -1;

        // count all and specific type instance in
        // all samples
        const std::array<std::string, 2> names {"all", type + ":" + kind + ":" + tag_string};
        for (auto const & name : names)
        {
            int i = 0;
            for (auto const &s : _impl->samples) {
                std::string key = name + ":" + s;

                // see if we already have a statistics counter for this kind of variant
                auto it = _impl->count_map.find(key);
                if (it == _impl->count_map.end()) {
                    it = _impl->count_map.emplace(key, VariantStatistics(*(_impl->fasta_to_use),
                                                                         _impl->count_homref)).first;
                }

                // determine the types seen in the variant
                int *types;
                int ntypes = 0;

                // count this variant
                it->second.add(_impl->hdr, v, i, &types, &ntypes);

                // aggregate the things we have seen across samples
                uint64_t vts_seen = 0;
                uint64_t lt = (uint64_t) -1;
                for (int j = 0; j < ntypes; ++j) {
                    vts_seen |= types[j] & 0xf;
                    vtypes.insert(types[j]);
                    if ((types[j] & 0xf0) > 0x10) {
                        lt = (uint64_t) (types[j] >> 4);
                    }
                }

                if (s == "TRUTH") {
                    hl_vt_truth = vts_seen == VT_NOCALL ? 2 : ((vts_seen == variant::VT_SNP
                                                                || vts_seen == (variant::VT_SNP | variant::VT_REF)) ? 0 : 1);
                    hl_lt_truth = lt;
                } else if (s == "QUERY") {
                    hl_vt_query = vts_seen == VT_NOCALL ? 2 : ((vts_seen == variant::VT_SNP
                                                                || vts_seen == (variant::VT_SNP | variant::VT_REF)) ? 0 : 1);
                    hl_lt_query = lt;
                }
                ++i;
            }
        }

        static const char *nvs[] = {"UNK", "SNP", "INDEL", "NOCALL"};
        if (hl_vt_truth < 4) {
            bcf_update_info_string(_impl->hdr, v, "T_VT", nvs[hl_vt_truth + 1]);
        }
        if (hl_vt_query < 4) {
            bcf_update_info_string(_impl->hdr, v, "Q_VT", nvs[hl_vt_query + 1]);
        }
        if (hl_lt_truth < 0x10) {
            bcf_update_info_string(_impl->hdr, v, "T_LT", CT_NAMES[hl_lt_truth]);
        }
        if (hl_lt_query < 0x10) {
            bcf_update_info_string(_impl->hdr, v, "Q_LT", CT_NAMES[hl_lt_query]);
        }
        if (_impl->output_vtc && !vtypes.empty()) {
            std::string s = "";
            int j = 0;
            for (int t : vtypes) {
                if (j++ > 0) { s += ","; }
                s += VariantStatistics::type2string(t);
            }
            bcf_update_info_string(_impl->hdr, v, "VTC", s.c_str());
        }
    }
}
