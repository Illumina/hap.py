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
 * Compares alleles in blocks
 *
 * \file BlockAlleleCompare.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <htslib/vcf.h>
#include "BlockAlleleCompare.hh"
#include "AlleleMatcher.hh"
#include "DistanceBasedMatcher.hh"

#include "helpers/BCFHelpers.hh"


namespace variant {
    struct BlockAlleleCompare::BlockAlleleCompareImpl
    {
        bcfhelpers::p_bcf_hdr hdr;
        FastaFile ref_fasta;
        std::string qq_name;
        int qq_hdr_id;
        bool qq_name_is_info;
        bool qq_name_is_qual;

        ComparisonMode comparison_mode = ALLELES;
        Json::Value comparison_parameters;

        std::list<bcfhelpers::p_bcf1> buffered;

        // track chrom / pos
        int current_rid = -1;
        int start_pos = -1;
    };

    /**
     * hdr must be destroyed externally, the reason it's not const is
     * that htslib doesn't do const mostly.
     *
     * @param hdr a bcf header
     * @param ref_fasta reference fasta file for trimming and counting
     */
    BlockAlleleCompare::BlockAlleleCompare(bcfhelpers::p_bcf_hdr hdr,
                                           FastaFile const & ref_fasta,
                                           std::string const & qq_name) : _impl(new BlockAlleleCompareImpl())
    {
        _impl->hdr = hdr;
        _impl->ref_fasta = ref_fasta;
        _impl->qq_name = qq_name;

        if(!qq_name.empty())
        {
            _impl->qq_hdr_id = bcf_hdr_id2int(hdr.get(), BCF_DT_ID, qq_name.c_str());
            if(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, _impl->qq_hdr_id))
            {
                _impl->qq_name_is_qual = false;
                _impl->qq_name_is_info = true;
            }
            else if(qq_name == "QUAL")
            {
                _impl->qq_name_is_qual = true;
                _impl->qq_name_is_info = false;
            }
            else
            {
                _impl->qq_name_is_info = false;
                _impl->qq_name_is_qual = false;
            }
        }
    }
    BlockAlleleCompare::~BlockAlleleCompare() {}

    BlockAlleleCompare::BlockAlleleCompare(BlockAlleleCompare && rhs)
    {
        _impl = std::move(rhs._impl);
    }
    BlockAlleleCompare & BlockAlleleCompare::operator=(BlockAlleleCompare && rhs)
    {
        _impl = std::move(rhs._impl);
        return *this;
    }

    /**
     * operators so we can sort blocks
     */
    bool BlockAlleleCompare::operator==(BlockAlleleCompare const & rhs) const
    {
        return _impl.get() == rhs._impl.get();
    }

    bool BlockAlleleCompare::operator<(BlockAlleleCompare const & rhs) const
    {
        return _impl->current_rid < rhs._impl->current_rid ||
               (_impl->current_rid == rhs._impl->current_rid && _impl->start_pos < rhs._impl->start_pos);
    }

    /**
     * Set the comparison mode
     */
    void BlockAlleleCompare::setComparisonMode(ComparisonMode mode)
    {
        _impl->comparison_mode = mode;
        _impl->comparison_parameters = Json::Value();
        switch(mode)
        {
            case ALLELES:
                break;
            case DISTANCE:
                _impl->comparison_parameters["max_distance"] = 30;
                break;
            default:
                break;
        }
    }

    /**
     * Set comparison parameters
     * @param params parameters for the comparison method
     */
    void BlockAlleleCompare::setComparisonParameters(Json::Value const & params)
    {
        _impl->comparison_parameters = params;
    }

    /**
     * Add a BCF record. Will duplicate the record and keep the copy
     * @param v bcf record
     */
    void BlockAlleleCompare::add(bcf1_t * v)
    {
        // make sure we don't change chromosome within a block
        if(_impl->current_rid < 0)
        {
            _impl->current_rid = v->rid;
        }
        else
        {
            assert(_impl->current_rid == v->rid);
        }
        _impl->buffered.emplace_back(bcfhelpers::pb(bcf_dup(v)));
        if(_impl->start_pos < 0 || v->pos < _impl->start_pos)
        {
            _impl->start_pos = v->pos;
        }
    }

    /**
     * Compare all buffered variants
     */
    void BlockAlleleCompare::run()
    {
        typedef std::list<std::pair<RefVar, bcf1_t*> > allele_list;
        allele_list truth_alleles;
        allele_list query_alleles;

        std::multimap<std::string, allele_list::iterator > truth_mapped;

        std::string current_chr;

        auto compare_and_update = [this,
                                   &current_chr,
                                   &truth_alleles,
                                   &query_alleles,
                                   &truth_mapped]()
        {
            std::unique_ptr<HapSetMatcher> p_matcher;

            switch(_impl->comparison_mode)
            {
                case ALLELES:
                    p_matcher.reset(new AlleleMatcher (_impl->ref_fasta.getFilename(),
                                                       current_chr));
                    break;
                case DISTANCE:
                    p_matcher.reset(new DistanceBasedMatcher (_impl->ref_fasta.getFilename(),
                                                              current_chr,
                                                              _impl->comparison_parameters["max_distance"].asInt()
                    ));
                    break;
                case ENUMERATE_DIPLOID:
                    p_matcher.reset(new HapSetMatcher(_impl->ref_fasta.getFilename(),
                                                      current_chr,
                                                      2,
                                                      _impl->comparison_parameters["max_enum"].asInt()
                                                      ));
                    break;
                default:
                    error("Invalid comparison mode: %i", _impl->comparison_mode);
            }

            std::unordered_map<size_t, bcf1_t*> comparison_mapping;
            std::unordered_map<size_t, int>     comparison_mapping_side;

            for(auto x = truth_alleles.begin(); x != truth_alleles.end(); ++x)
            {
                const size_t v_id = p_matcher->addLeft(x->first, 1);
                comparison_mapping[v_id] = x->second;
                comparison_mapping_side[v_id] = 0;
            }

            for(auto x = query_alleles.begin(); x != query_alleles.end(); ++x)
            {
                const size_t v_id = p_matcher->addRight(x->first, 1);
                comparison_mapping[v_id] = x->second;
                comparison_mapping_side[v_id] = 1;
            }

            HapAssignment assignment;
            p_matcher->optimize(assignment);

            struct Updates {
                std::vector<std::string> bd = {".", "."};
            };

            std::map<bcf1_t*, Updates> update_list;
            auto u = [&update_list](bcf1_t* v) -> std::map<bcf1_t*, Updates>::iterator
            {
                auto it = update_list.find(v);
                if(it == update_list.end())
                {
                    it = update_list.emplace(std::make_pair(v, Updates())).first;
                }
                return it;
            };

            for(size_t v_id = 0; v_id < assignment.variant_assignments.size(); ++v_id)
            {
                if(assignment.variant_assignments[v_id] != 0)
                {
                    u(comparison_mapping[v_id])->second.bd[comparison_mapping_side[v_id]] = "TP";
                }
                else if(comparison_mapping_side[v_id] == 0)
                {
                    u(comparison_mapping[v_id])->second.bd[0] = "FN";
                }
                else
                {
                    u(comparison_mapping[v_id])->second.bd[1] = "FP";
                }
            }

            for(auto & update : update_list)
            {
                bcfhelpers::setFormatStrings(_impl->hdr.get(), update.first, "BD", update.second.bd);
            }

            query_alleles.clear();
            truth_alleles.clear();
        };

        int64_t current_block_end = -1;

        for(auto & v : _impl->buffered)
        {
            bcf_unpack(v.get(), BCF_UN_ALL);

            bool phased = false;
            int ngt = 0;
            int gt[MAX_GT];
            const std::string vchr = bcfhelpers::getChrom(_impl->hdr.get(), v.get());

            // make sure all variants in block are on same chr
            if(!current_chr.empty() && current_chr != vchr)
            {
                compare_and_update();
                current_block_end = -1;
            }
            current_chr = vchr;

            auto makeRefVar = [this, &vchr, &v](int allele) -> RefVar
            {
                assert(allele < v->n_allele);
                RefVar rv;
                rv.start = v->pos;
                rv.alt = v->d.allele[allele];
                if(rv.alt == "<DEL>")
                {
                    rv.alt = "";
                }
                else if(rv.alt[0] == '<')
                {
                    error("unsupported symbolic ALT at %s:%i : %S", vchr.c_str(), v->pos, rv.alt.c_str());
                }
                rv.end = rv.start + rv.alt.size() - 1;
                return rv;
            };

            bcfhelpers::getGT(_impl->hdr.get(), v.get(), 0, gt, ngt, phased);
            for(int igt = 0; igt < ngt; ++igt)
            {
                if(gt[igt] > 0)
                {
                    truth_alleles.push_back(std::make_pair(makeRefVar(gt[igt]), v.get()));
                }
            }

            bcfhelpers::getGT(_impl->hdr.get(), v.get(), 1, gt, ngt, phased);
            for(int igt = 0; igt < ngt; ++igt)
            {
                if(gt[igt] > 0)
                {
                    query_alleles.push_back(std::make_pair(makeRefVar(gt[igt]), v.get()));
                }
            }

            std::vector<float> QQ = {std::numeric_limits<float>::quiet_NaN(),
                                     std::numeric_limits<float>::quiet_NaN()};
            if(_impl->qq_name_is_info)
            {
                QQ[0] = QQ[1] = bcfhelpers::getInfoFloat(_impl->hdr.get(), v.get(), _impl->qq_name.c_str());
            }
            else if(_impl->qq_name_is_qual)
            {
                QQ[0] = QQ[1] = v->qual;
            }
            else
            {
                QQ[0] = bcfhelpers::getFormatFloat(_impl->hdr.get(), v.get(), _impl->qq_name.c_str(), 0);
                QQ[1] = bcfhelpers::getFormatFloat(_impl->hdr.get(), v.get(), _impl->qq_name.c_str(), 1);
            }
            bcfhelpers::setFormatFloats(_impl->hdr.get(), v.get(), "QQ", QQ);

            if(  (current_block_end >= 0 && v->pos > current_block_end + 100000)
               || truth_alleles.size() + query_alleles.size() > 1024)
            {
                compare_and_update();
                current_block_end = -1;
            }
            int64_t vstart, vend;
            bcfhelpers::getLocation(_impl->hdr.get(), v.get(), vstart, vend);
            current_block_end = std::max(vstart, current_block_end);
            current_block_end = std::max(vend, current_block_end);
        }

        if(!current_chr.empty() && !(truth_alleles.empty() && query_alleles.empty()))
        {
            compare_and_update();
        }
    }

    /**
     * Output all buffered variants into a file
     * @param output HTS file to write to
     */
    void BlockAlleleCompare::output(htsFile * output)
    {
        for(auto & v : _impl->buffered)
        {
            bcf_write(output, _impl->hdr.get(), v.get());
        }
    }

    /**
     * BlockAlleleCompare::run doesn't write a header such that
     * results can be concatenated. This function updates the header once in the beginning
     * @param hdr header to update
     */
    void BlockAlleleCompare::updateHeader(bcf_hdr_t * hdr)
    {
        if(bcf_hdr_nsamples(hdr) != 2)
        {
            error("Input file must have exactly two samples. We have %i", bcf_hdr_nsamples(hdr));
        }

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
        bcf_hdr_append(hdr, "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">");

        bcf_hdr_sync(hdr);
    }
}


