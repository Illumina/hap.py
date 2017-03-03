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
        auto p1 = params.find("QQ:");
        if(p1 != std::string::npos)
        {
            std::string qqstr = params.substr(p1 + 3);
            auto p2 = qqstr.find(";");
            if(p2 != std::string::npos)
            {
                qqstr = qqstr.substr(0, p2);
            }
            roc_field = qqstr;
        }

        if(!roc_field.empty())
        {
            roc_hdr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, roc_field.c_str());
            if(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, roc_hdr_id))
            {
                roc_field_is_qual = false;
                roc_field_is_info = true;
            }
            else if(roc_field == "QUAL")
            {
                roc_field_is_qual = true;
                roc_field_is_info = false;
            }
            else
            {
                roc_field_is_info = false;
                roc_field_is_qual = false;
            }
        }

        p1 = params.find("clean_info");
        clean_info = p1 != std::string::npos;
    }

    void XCMPQuantify::updateHeader(bcf_hdr_t * hdr)
    {
        // check if we have all the required fields present
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "BS", NULL))
        {
            bcf_hdr_append(hdr, "##INFO=<ID=BS,Number=.,Type=Integer,Description=\"Benchmarking superlocus ID for these variants.\">");
        }
        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "GT", NULL))
        {
            bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
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
            bcf_hdr_append(hdr, "##INFO=<ID=XCMP,Number=.,Type=String,Description=\"XCMP extra information.\">");
        }
        bcf_hdr_sync(hdr);
    }

    void XCMPQuantify::countVariants(bcf1_t * v)
    {
        bcf_unpack(v, BCF_UN_ALL);

        const std::string tag_string = bcfhelpers::getInfoString(_impl->hdr, v, "Regions", "");
        std::string type = bcfhelpers::getInfoString(_impl->hdr, v, "type");
        std::string kind = bcfhelpers::getInfoString(_impl->hdr, v, "kind");
        std::string ctype = bcfhelpers::getInfoString(_impl->hdr, v, "ctype");
        std::string gtt1 = ".";
        std::string gtt2 = ".";

        if(_impl->output_vtc)
        {
            gtt1 = bcfhelpers::getInfoString(_impl->hdr, v, "gtt1");
            gtt2 = bcfhelpers::getInfoString(_impl->hdr, v, "gtt2");
        }
        const bool hapmatch = bcfhelpers::getInfoFlag(_impl->hdr, v, "HapMatch");
        const bool fail = bcfhelpers::getInfoFlag(_impl->hdr, v, "IMPORT_FAIL");
        const bool q_filtered = bcfhelpers::getInfoFlag(_impl->hdr, v, "Q_FILTERED");
        float QQ = std::numeric_limits<float>::quiet_NaN();
        if(roc_field_is_info)
        {
            QQ = bcfhelpers::getInfoFloat(_impl->hdr, v, roc_field.c_str());
        }
        else if(roc_field_is_qual)
        {
            QQ = v->qual;
        }

        /** clear all info fields except for GA4GH fields and END */
        if(clean_info)
        {
            for (int i = 0; i < v->n_info; ++i) {
                const char * key = bcf_hdr_int2id(_impl->hdr, BCF_DT_ID, v->d.info[i].key);
                if(    strcmp(key, "END") != 0
                    && strcmp(key, "VTC") != 0
                    && strcmp(key, "Regions") != 0
                    && strcmp(key, "BS") != 0
                    && strcmp(key, "XCMP") != 0
                    && strcmp(key, "IMPORT_FAIL") != 0
                )
                {
                    const int ntype = v->d.info[i].type;
                    bcf_update_info(_impl->hdr, v, key, NULL, 0, ntype);
                }
            }
        }

        if (hapmatch && type != "TP" && !q_filtered) {
            kind = "hapmatch__" + type + "__" + kind;
            type = "TP";
        }

        if (_impl->count_unk && tag_string.find("CONF") == std::string::npos) {
            type = "UNK";
        }

        if(fail)
        {
            type = "N";
            kind = "error";
        }

        if(_impl->output_vtc)
        {
            bcf_update_info_string(_impl->hdr, v, "XCMP",
                                   (type + ":" + kind + ":" + gtt1 + ":" + gtt2 + ":" + ctype).c_str());
        }

        if(type == "TP")
        {
            kind = "gm";
        }
        else if(kind == "gtmismatch")
        {
            kind = "am";
        }
        else if(kind == "almismatch" || ctype == "hap:mismatch")
        {
            kind = "lm";
        }
        else
        {
            kind = ".";
        }

        std::set<int> vtypes;
        std::vector<std::string> bds;
        std::vector<std::string> bks;
        std::vector<std::string> bis;
        std::vector<std::string> vts;
        std::vector<std::string> lts;
        std::vector<float> qqs;

        // count all and specific type instance in
        // all samples
        int i = 0;
        for (auto const &s : _impl->samples) {
            // figure out if this is a no-call first
            int gt[MAX_GT];
            int ngt = 0;
            bool phased = false;
            bcfhelpers::getGT(_impl->hdr, v, i, gt, ngt, phased);
            const bool isNoCall = ngt == 0 || std::all_of(&gt[0], &gt[ngt], [](int x) { return x < 0; });
            const bool isHomref = ngt > 0 && std::all_of(&gt[0], &gt[ngt], [](int x) { return x == 0; });

            // count as "all"
            std::string key = "all:" + s;

            // see if we already have a statistics counter for this kind of variant
            auto it = _impl->count_map.find(key);
            if (it == _impl->count_map.end()) {
                it = _impl->count_map.emplace(key, VariantStatistics(*(_impl->fasta_to_use),
                                                                     _impl->count_homref)).first;
            }
            int *types;
            int ntypes = 0;
            std::set<std::string> * p_extracounts = NULL;
            it->second.add(_impl->hdr, v, i, &types, &ntypes, &p_extracounts);

            // aggregate the things we have seen across samples
            uint64_t vts_seen = 0;
            uint64_t lt = (uint64_t) -1;
            for (int j = 0; j < ntypes; ++j)
            {
                vts_seen |= types[j] & 0xf;
                vtypes.insert(types[j]);
                if ((types[j] & 0xf0) > 0x10)
                {
                    lt = (uint64_t) (types[j] >> 4);
                }
            }

            if(p_extracounts)
            {
                bis.push_back(VariantStatistics::extraCountsToBI(*p_extracounts));
            }
            else
            {
                bis.push_back(".");
            }

            /** vt = 0 -> UNK
             *  vt = 1 -> SNP
             *  vt = 2 -> INDEL
             *  vt = 3 -> NOCALL
             *  vt = 4 -> HOMREF
             */
            static const char *nvs[] = {"UNK", "SNP", "INDEL", "NOCALL", "HOMREF"};
            uint64_t vt = 0;

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

            // determine the types seen in the variant
            if(fail || isNoCall || (i == 1 && q_filtered))
            {
                key = ".:.:" + tag_string + ":" + s;

                it = _impl->count_map.find(key);
                if (it == _impl->count_map.end()) {
                    it = _impl->count_map.emplace(key, VariantStatistics(*(_impl->fasta_to_use),
                                                                         _impl->count_homref)).first;
                }
                // count this variant
                it->second.add(_impl->hdr, v, i);

                if(fail)
                {
                    bds.push_back("N");
                }
                else
                {
                    bds.push_back(".");
                }
                bks.push_back(".");
                vts.push_back("NOCALL");
                if(fail)
                {
                    lts.push_back("fail");
                }
                else if(q_filtered)
                {
                    lts.push_back("filtered");
                }
                else
                {
                    lts.push_back("nocall");
                }
                qqs.push_back(0);
            }
            else
            {
                std::string this_type = type;
                /** don't annotate no-calls */
                if(vts_seen == VT_NOCALL)
                {
                    this_type = "N";
                    bds.push_back(".");
                    bks.push_back(".");
                }
                else
                {
                    if(i == 0 && type == "FP")
                    {
                        // GT/allele mismatches become FNs
                        bds.push_back("FN");
                        this_type = "FN";
                    }
                    else
                    {
                        bds.push_back(type);
                    }
                    bks.push_back(kind);
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

                // determine QQ
                if(roc_field_is_info || roc_field_is_qual)
                {
                    qqs.push_back(QQ);
                }
                else
                {
                    qqs.push_back(bcfhelpers::getFormatFloat(_impl->hdr, v, roc_field.c_str(), i));
                }

                key = this_type + ":" + kind + ":" + tag_string + ":" + s;

                it = _impl->count_map.find(key);
                if (it == _impl->count_map.end())
                {
                    it = _impl->count_map.emplace(key, VariantStatistics(*(_impl->fasta_to_use),
                                                                         _impl->count_homref)).first;
                }

                // count this variant
                it->second.add(_impl->hdr, v, i);
            }
            ++i;
        }

        bcfhelpers::setFormatStrings(_impl->hdr, v, "BD", bds);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BK", bks);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BI", bis);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BVT", vts);
        bcfhelpers::setFormatStrings(_impl->hdr, v, "BLT", lts);
        bcfhelpers::setFormatFloats(_impl->hdr, v, "QQ", qqs);

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
