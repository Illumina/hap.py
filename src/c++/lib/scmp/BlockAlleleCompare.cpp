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

        std::list<bcfhelpers::p_bcf1> buffered;
    };

    /**
     * hdr must be destroyed externally, the reason it's not const is
     * that htslib doesn't do const mostly.
     *
     * @param hdr a bcf header
     * @param ref_fasta reference fasta file for trimming and counting
     */
    BlockAlleleCompare::BlockAlleleCompare(bcf_hdr_t * hdr,
                                           FastaFile const & ref_fasta,
                                           std::string const & qq_name) : _impl(new BlockAlleleCompareImpl())
    {
        _impl->hdr = bcfhelpers::ph(bcf_hdr_dup(hdr));
        _impl->ref_fasta = ref_fasta;
        _impl->qq_name = qq_name;

        if(!qq_name.empty())
        {
            _impl->qq_hdr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, qq_name.c_str());
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
     * Add a BCF record. Will duplicate the record and keep the copy
     * @param v bcf record
     */
    void BlockAlleleCompare::add(bcf1_t * v)
    {
        _impl->buffered.emplace_back(bcfhelpers::pb(bcf_dup(v)));
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
        std::list< std::pair<bcf1_t*, bcf1_t*> > direct_matches;

        std::list< allele_list::iterator > truth_unmatched;
        std::list< allele_list::iterator > query_unmatched;

        auto compare_and_update = [&_impl,
                                   &truth_alleles,
                                   &query_alleles,
                                   &direct_matches,
                                   &truth_unmatched,
                                   &query_unmatched](bool finish)
        {

            for(auto x = truth_alleles.begin(); x != truth_alleles.end(); ++x)
            {
                truth_mapped.insert(std::make_pair(x->first.repr(), x));
            }

            for(auto x = query_alleles.begin(); x != query_alleles.end(); ++x)
            {
                auto it = truth_mapped.find(x->first.repr());
                if(it != truth_mapped.end())
                {
                    direct_matches.emplace_back(std::make_pair(it->second->second, x->second));
                    truth_mapped.erase(it);
                }
                else
                {
                    query_unmatched.push_back(x);
                }
            }
            for(auto & x : truth_mapped)
            {
                truth_unmatched.push_back(x.second);
            }

            if(finish || (truth_unmatched.empty() && query_unmatched.empty()))
            {
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

                for(auto & match : direct_matches)
                {
                    u(match.first)->second.bd[0] = "TP";
                    u(match.second)->second.bd[1] = "TP";
                }

                for(auto & fn : truth_unmatched)
                {
                    u(fn->second)->second.bd[0] = "FN";
                }

                for(auto & fp : query_unmatched)
                {
                    u(fp->second)->second.bd[1] = "FP";
                }

                for(auto & update : update_list)
                {
                    bcfhelpers::setFormatStrings(_impl->hdr.get(), update.first, "BD", update.second.bd);
                }

                query_alleles.clear();
                truth_alleles.clear();
                truth_unmatched.clear();
                query_unmatched.clear();
                direct_matches.clear();
            }
        };


        for(auto & v : _impl->buffered)
        {
            bcf_unpack(v.get(), BCF_UN_ALL);

            bool phased = false;
            int ngt = 0;
            int gt[MAX_GT];
            const std::string vchr = bcfhelpers::getChrom(_impl->hdr.get(), v.get());

            auto makeRefVar = [&_impl, &vchr, &v](int allele) -> RefVar
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
                trimLeft(_impl->ref_fasta, vchr.c_str(), rv, false);
                trimRight(_impl->ref_fasta, vchr.c_str(), rv, false);
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

            to_update.push_back(v.get());
            compare_and_update(false);
        }

        compare_and_update(true);
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


