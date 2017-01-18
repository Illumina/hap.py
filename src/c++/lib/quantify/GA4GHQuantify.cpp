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
 * GA4GH GT match quantification
 *
 * See
 * https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md
 *
 * \file GA4GHQuantify.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "GA4GHQuantify.hh"

/* #define DEBUG_GA4GH_QUANTIFY */

namespace variant
{

    GA4GHQuantify::GA4GHQuantify(bcf_hdr_t *hdr,
                               FastaFile const &ref_fasta,
                               const std::string & params) :
        BlockQuantify(hdr, ref_fasta, params)
    {

    }

    void GA4GHQuantify::updateHeader(bcf_hdr_t *hdr)
    {
        // check if we have all the required fields present
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "BS", NULL))
        {
            bcf_hdr_append(hdr, "##INFO=<ID=BS,Number=.,Type=Integer,Description=\"Benchmarking superlocus ID for these variants.\">");
        }
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "GT", NULL))
        {
            bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        }
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "BD", NULL))
        {
            bcf_hdr_append(hdr, "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
        }
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "BK", NULL))
        {
            bcf_hdr_append(hdr,
                           "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
        }
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "BI", NULL))
        {
            bcf_hdr_append(hdr, "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">");
        }
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "QQ", NULL))
        {
            bcf_hdr_append(hdr, "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">");
        }
        bcf_hdr_append(hdr, "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">");
        if(_impl->output_vtc)
        {
            bcf_hdr_append(hdr, "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">");
        }
        bcf_hdr_sync(hdr);
    }

    void GA4GHQuantify::countVariants(bcf1_t * v)
    {
        bcf_unpack(v, BCF_UN_ALL);
        std::string tag_string = bcfhelpers::getInfoString(_impl->hdr, v, "Regions", "");
        std::set<int> vtypes;
        std::vector<std::string> bds;
        std::vector<std::string> bis;
        std::vector<std::string> vts;
        std::vector<std::string> lts;

        int si = 0;
        for (auto const &s : _impl->samples) {
            std::string type = bcfhelpers::getFormatString(_impl->hdr, v, "BD", si);
            std::string kind = bcfhelpers::getFormatString(_impl->hdr, v, "BK", si);

            if (_impl->count_unk && tag_string.find("CONF") == std::string::npos) {
                type = "UNK";
            }

            bds.push_back(type);

            // count "all"
            std::string key = "all:" + s;
            auto it = _impl->count_map.find(key);
            if (it == _impl->count_map.end()) {
                it = _impl->count_map.emplace(key, VariantStatistics(*(_impl->fasta_to_use),
                                                                     _impl->count_homref)).first;
            }
            it->second.add(_impl->hdr, v, si);

            // count all and specific type instance in
            // all samples
            key = type + ":" + kind + ":" + tag_string + ":" + s;
#ifdef DEBUG_GA4GH_QUANTIFY
            std::cerr << key << "\n";
#endif

            it = _impl->count_map.find(key);
            if (it == _impl->count_map.end()) {
                it = _impl->count_map.emplace(key, VariantStatistics(*(_impl->fasta_to_use),
                                                                     _impl->count_homref)).first;
            }

            // determine the types seen in the variant
            int *types;
            int ntypes = 0;
            std::set<std::string> * p_extracounts = NULL;

            // count this variant
            it->second.add(_impl->hdr, v, si, &types, &ntypes, &p_extracounts);

            if(p_extracounts)
            {
                bis.push_back(VariantStatistics::extraCountsToBI(*p_extracounts));
            }
            else
            {
                bis.push_back(".");
            }

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

            /** vt = 0 -> UNK
             *  vt = 1 -> SNP
             *  vt = 2 -> INDEL
             *  vt = 3 -> NOCALL
             *  vt = 4 -> HOMREF
             */
            static const char *nvs[] = {"UNK", "SNP", "INDEL", "NOCALL", "HOMREF"};
            uint64_t vt = 0;
            int gt[MAX_GT];
            int ngt = 0;
            bool phased = false;
            bcfhelpers::getGT(_impl->hdr, v, si, gt, ngt, phased);
            const bool isNoCall = ngt == 0 || std::all_of(&gt[0], &gt[ngt], [](int x) { return x < 0; });
            const bool isHomref = ngt > 0 && std::all_of(&gt[0], &gt[ngt], [](int x) { return x == 0; });

            if(vts_seen == variant::VT_SNP || vts_seen == (variant::VT_SNP | variant::VT_REF))
            {
                vt = 1;
            }
            else if(vts_seen & VT_INS || vts_seen & VT_DEL || vts_seen & VT_SNP)
            {
                vt = 2;
            }
            if(isNoCall)
            {
                vt = 3;
            }

            if(isHomref)
            {
                vt = 4;
            }
            vts.push_back(nvs[vt]);
            if(lt < 0x0f)
            {
                lts.push_back(CT_NAMES[lt]);
            }
            else
            {
                lts.push_back("error");
            }
            ++si;
        }

        bcfhelpers::setFormatStrings(_impl->hdr, v, "BD", bds);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BI", bis);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BVT", vts);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BLT", lts);
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
