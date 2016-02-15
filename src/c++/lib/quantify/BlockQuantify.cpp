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

#include "BlockQuantify.hh"


//#define DEBUG_BLOCKQUANTIFY

namespace variant {

    struct BlockQuantify::BlockQuantifyImpl {
        ~BlockQuantifyImpl()
        {
            for(auto x : variants)
            {
                bcf_destroy(x);
            }
        }

        bcf_hdr_t * hdr;
        FastaFile const & ref_fasta;
        std::unique_ptr<FastaFile> fasta_to_use;
        bool output_vtc;
        bool count_homref;
        bool count_unk;

        typedef std::list<std::string> samplenames_t;
        typedef std::map<std::string, VariantStatistics> count_map_t;
        typedef std::list<bcf1_t *> variantlist_t;

        count_map_t count_map;
        variantlist_t variants;
        samplenames_t samples;
    };

    BlockQuantify::BlockQuantify(bcf_hdr_t * hdr,
                                 FastaFile const & ref_fasta,
                                 bool output_vtc,
                                 bool count_homref,
                                 bool count_unk) :
        _impl(std::unique_ptr<BlockQuantifyImpl>(new BlockQuantifyImpl{hdr,
                                                    ref_fasta,
                                                    nullptr,
                                                    output_vtc,
                                                    count_homref,
                                                    count_unk,
                                                    BlockQuantifyImpl::count_map_t(),
                                                    BlockQuantifyImpl::variantlist_t(),
                                                    bcfhelpers::getSampleNames(hdr)}))
    {
    }

    BlockQuantify::~BlockQuantify() {}

    BlockQuantify::BlockQuantify(BlockQuantify && rhs) : _impl(std::move(rhs._impl)) { }

    BlockQuantify & BlockQuantify::operator=(BlockQuantify && rhs)
    {
        _impl = std::move(rhs._impl);
        return *this;
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
#endif
        for(auto & v : _impl->variants)
        {
#ifdef DEBUG_BLOCKQUANTIFY
            lastpos = v->pos;
#endif
            count_variants(v);
        }
#ifdef DEBUG_BLOCKQUANTIFY
        std::cerr << "finished block " << lastpos << " - " << _impl->variants.size() << " records on thread " << std::this_thread::get_id() << "\n";
#endif
        _impl->fasta_to_use.reset(nullptr);
    }

    // result output
    std::map<std::string, VariantStatistics> const & BlockQuantify::getCounts() const
    {
        return _impl->count_map;
    }

    std::list<bcf1_t*> const & BlockQuantify::getVariants()
    {
        return _impl->variants;
    }

    void BlockQuantify::count_variants(bcf1_t * v) {
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
