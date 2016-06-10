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
 *  \brief Variant I/O implementation
 *
 * \file VariantWriter.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "VariantImpl.hh"

#include "helpers/StringUtil.hh"

#include <sstream>
#include <cstring>
#include <unordered_set>
#include <htslib/vcf.h>

// #define DEBUG_VARIANT_GTS

namespace variant
{

    VariantWriter::VariantWriter(const char * filename, const char * reference)
    {
        _impl = new VariantWriterImpl(filename, reference);
    }

    VariantWriter::VariantWriter(VariantWriter const & rhs)
    {
        _impl = new VariantWriterImpl(rhs._impl->filename.c_str(), rhs._impl->referencename.c_str());
        for(std::string const & s : rhs._impl->samples)
        {
            addSample(s.c_str());
        }
        _impl->header_lines = rhs._impl->header_lines;
    }
    VariantWriter::~VariantWriter()
    {
        if(!_impl->header_done)
        {
            _impl->writeHeader();
        }
        delete _impl;
    }

    VariantWriter const & VariantWriter::operator=(VariantWriter const & rhs)
    {
        if(&rhs == this)
        {
            return *this;
        }
        delete _impl;
        _impl = new VariantWriterImpl(rhs._impl->filename.c_str(), rhs._impl->referencename.c_str());
        for(std::string const & s : rhs._impl->samples)
        {
            addSample(s.c_str());
        }
        _impl->header_lines = rhs._impl->header_lines;
        return *this;
    }

    /** write format fields apart from GT */
    void VariantWriter::setWriteFormats(bool write_fmts)
    {
        _impl->write_formats = write_fmts;
    }

    bool VariantWriter::getWriteFormats() const
    {
        return _impl->write_formats;
    }

    void VariantWriter::addHeader(const char * headerline)
    {
        _impl->header_lines.push_back(headerline);
    }

    void VariantWriter::addHeader(VariantReader & vr, bool drop)
    {
        for (int i = 0; i < vr._impl->files->nreaders; ++i)
        {
            std::vector<std::string> result;
            int len = 0;
            char * hdr_text = bcf_hdr_fmt_text(vr._impl->files->readers[i].header, 0, &len);
            if(!hdr_text)
            {
                continue;
            }
            stringutil::split(std::string(hdr_text, len), result, "\n");
            free(hdr_text);
            for(std::string hl : result) {
                if(drop && !(
                    (hl.find("##contig") == 0) ||
                    (hl.find("##FILTER") == 0) ||
                    (hl.find("##ref") == 0) )
                ) { continue; }
                _impl->header_lines.push_back(hl);
            }
        }
    }

    int VariantWriter::addSample(const char * sname)
    {
        _impl->samples.push_back(std::string(sname));

        bcf_hdr_add_sample(_impl->hdr, sname);

        return (int)_impl->samples.size()-1;
    }

    void VariantWriterImpl::writeHeader()
    {
        // merge header permissively
        std::map<std::string, std::string> merged_header_items;

        std::sort(header_lines.begin(), header_lines.end());

        for (std::string const & l : header_lines)
        {
            int len = 0;
            if(l.size() < 2 || l[0] != '#' || l[1] != '#')
            {
                continue;
            }
            bcf_hrec_t * hrec = bcf_hdr_parse_line(hdr, l.c_str(), &len);
            if(!hrec)
            {
                std::cerr << "[W] Failed to parse VCF header line: " << l << "\n";
                continue;
            }
            int idx = bcf_hrec_find_key(hrec, "ID");
            if ( idx < 0 )
            {
                // no ID => skip
                bcf_hrec_destroy(hrec);
                continue;
            }

            std::ostringstream oss;
            oss << hrec->key << ":";
            if(hrec->value)
            {
                oss << hrec->value << ":";
            }
            else
            {
                oss << ":";
            }
            oss << hrec->vals[idx];

            if(merged_header_items.find(oss.str()) == merged_header_items.end())
            {
                merged_header_items[oss.str()] = l;
                if(hrec->type == BCF_HL_INFO || l.substr(0, 6) == "##INFO") {
                    // HAP-79 check if info entries have a description
                    int idx = bcf_hrec_find_key(hrec, "Description");
                    if (idx < 0)
                    {
                        bcf_hrec_add_key(hrec, "Description", 11);
                        idx = bcf_hrec_find_key(hrec, "Description");
                        bcf_hrec_set_val(hrec, idx, "none", 4, 1);
                    }
                }
                bcf_hdr_add_hrec(hdr, hrec);
            }
            else
            {
                // skip duplicate
                bcf_hrec_destroy(hrec);
            }
        }

        bcf_hdr_add_sample(hdr, NULL);
        bcf_hdr_set_version(hdr, "VCFv4.1");
        bcf_hdr_write(fp, hdr);
        header_done = true;
    }

    void VariantWriter::put(Variants const & var)
    {
        if(!_impl->header_done)
        {
            _impl->writeHeader();
        }

        bcf_hdr_t * hdr = _impl->hdr;
        bcf1_t *rec = _impl->rec;
        bcf_clear(rec);

        rec->rid = bcf_hdr_name2id(hdr, var.chr.c_str());
        rec->pos = var.pos;

        /** make alleles */
        if(var.variation.size() > 0)
        {
            std::vector<std::string> alleles;
            toAlleles(_impl->reference, var.chr.c_str(), var.variation, alleles);

            std::ostringstream ss;
            for(size_t j = 0; j < alleles.size(); ++j)
            {
                if(j > 0)
                {
                    ss << ",";
                }
                ss << alleles[j];
            }
            bcf_update_alleles_str(hdr, rec, ss.str().c_str());
        }
        else    // homref?
        {
            std::string refnt = _impl->reference.query(var.chr.c_str(), rec->pos, rec->pos);
            bcf_update_alleles_str(hdr, rec, refnt.c_str());
            int varend = (int) (var.pos + var.len);
            bcf_update_info_int32(hdr, rec, "END", &varend, 1);
        }

        // .. QUAL => max over all calls
        rec->qual = 0;
        for(Call const & c : var.calls)
        {
            if(!c.bcf_rec)
            {
                continue;
            }
            const float cqual = c.bcf_rec->qual;
            if(!std::isnan(cqual))
            {
                rec->qual = std::max(cqual, rec->qual);
            }
        }

        int fmap[MAX_FILTER];
        int fcount = 0;
        std::set<int> added;
        memset(fmap, -1, sizeof(int)*MAX_FILTER);

        // combine all filters
        for(Call const & call : var.calls)
        {
            for(uint64_t c = 0; c < call.nfilter; ++c)
            {
                if(fcount >= MAX_FILTER)
                {
                    error("Record has too many filters at %s:%i", var.chr.c_str(), rec->pos);
                }
                int k = bcf_hdr_id2int(hdr, BCF_DT_ID, call.filter[c].c_str());
                if(k < 0)
                {
                    error("Filter '%s' is not in the VCF header. This will break VCF writing, at %s:%i",
                            call.filter[c].c_str(), var.chr.c_str(), rec->pos);
                }

                if(!added.count(k))
                {
                    // don't add PASS
                    if(call.filter[c] != "PASS" && call.filter[c] != ".")
                    {
                        fmap[fcount++] = k;
                        added.insert(k);
                    }
                }
            }
        }
        bcf_update_filter(hdr, rec, fmap, fcount);

        // write GTs
        int32_t *tmp = (int32_t*)malloc(bcf_hdr_nsamples(hdr)*MAX_GT*sizeof(int32_t));
        memset(tmp, 0, sizeof(int32_t)*bcf_hdr_nsamples(hdr)*MAX_GT);

        for(size_t g = 0; g < var.calls.size(); ++g)
        {
            Call const & c = var.calls[g];

            // optionally write other formats
            for (size_t j = 0; j < var.calls[g].ngt; ++j)
            {
                if(c.phased)
                {
                    tmp[j + MAX_GT*g] = bcf_gt_phased(c.gt[j]);
                }
                else
                {
                    tmp[j + MAX_GT*g] = bcf_gt_unphased(c.gt[j]);
                }
            }
        }

        if(bcf_update_genotypes(hdr, rec, tmp, bcf_hdr_nsamples(hdr)*MAX_GT) < 0)
        {
            error("Cannot update bcf genotypes at %s:%i", var.chr.c_str(), rec->pos);
        }

        free(tmp);

        // write info
        {
            // prevent keeping AN and AC since these might change
            std::set<std::string> infos_written = {"AN", "AC"};
            std::map<std::string, std::vector<int> > int_infos;
            std::map<std::string, std::vector<float> > float_infos;
            std::map<std::string, std::string > string_infos;
            for(auto & c : var.calls)
            {
                bcf_hdr_t * bhdr = c.bcf_hdr.get();
                bcf1_t * bline = c.bcf_rec.get();
                bcf_unpack(bline, BCF_UN_INFO);

                for(int ni = 0; ni < bline->n_info; ++ni)
                {
                    bcf_info_t * inf = &bline->d.info[ni];
                    const char * id = bcf_hdr_int2id(bhdr, BCF_DT_ID, inf->key);
                    if(infos_written.count(id))
                    {
                        continue;
                    }
                    switch(inf->type)
                    {
                        case BCF_BT_INT8:
                        case BCF_BT_INT16:
                        case BCF_BT_INT32:
                            int_infos.emplace(std::make_pair(
                                std::string(id),
                                bcfhelpers::getInfoInts(bhdr, bline, id)));
                            infos_written.insert(id);
                            break;
                        case BCF_BT_FLOAT:
                            float_infos.emplace(std::make_pair(
                                std::string(id),
                                bcfhelpers::getInfoFloats(bhdr, bline, id)));
                            infos_written.insert(id);
                            break;
                        case BCF_BT_CHAR:
                            string_infos.emplace(std::make_pair(
                                std::string(id),
                                bcfhelpers::getInfoString(bhdr, bline, id)));
                            infos_written.insert(id);
                            break;
                        default:
                            // default: ignore
                            break;
                    }
                }
            }

            for(auto const & x : int_infos)
            {
                auto p = std::unique_ptr<int[]>(new int[x.second.size()]);
                for(size_t j = 0; j < x.second.size(); ++j)
                {
                    p.get()[j] = x.second[j];
                }
                bcf_update_info_int32(hdr, rec, x.first.c_str(), p.get(), (int)x.second.size());
            }
            for(auto const & x : float_infos)
            {
                auto p = std::unique_ptr<float[]>(new float[x.second.size()]);
                for(size_t j = 0; j < x.second.size(); ++j)
                {
                    p.get()[j] = x.second[j];
                }
                bcf_update_info_float(hdr, rec, x.first.c_str(), p.get(), (int)x.second.size());
            }
            for(auto const & x : string_infos)
            {
                bcf_update_info_string(hdr, rec, x.first.c_str(), x.second.c_str());
            }
        }

        if (_impl->write_formats)
        {
            auto tmp_dp = std::unique_ptr<int32_t[]>(new int32_t [bcf_hdr_nsamples(hdr)]);
            auto tmp_ad = std::unique_ptr<int32_t[]>(new int32_t [bcf_hdr_nsamples(hdr) * (var.variation.size() + 1)]);
            auto tmp_ado = std::unique_ptr<int32_t[]>(new int32_t [bcf_hdr_nsamples(hdr)]);

            std::map<std::string, std::vector<int> > int_fmts;
            std::map<std::string, std::vector<float> > float_fmts;
            std::map<std::string, std::vector<std::string> > string_fmts;

            memset(tmp_ad.get(), 0, sizeof(int32_t)*bcf_hdr_nsamples(hdr)*(var.variation.size() + 1));
            memset(tmp_ado.get(), 0, sizeof(int32_t)*bcf_hdr_nsamples(hdr));
            memset(tmp_dp.get(), 0, sizeof(int32_t)*bcf_hdr_nsamples(hdr));

            for(size_t g = 0; g < var.calls.size(); ++g)
            {
                if(int(g) >= bcf_hdr_nsamples(hdr))
                {
                    std::cerr << "[W] Variant record contains too many calls for writer." << "\n";
                    break;
                }
                Call const & c = var.calls[g];

                tmp_dp[g] = c.dp < 0 ? bcf_int32_missing : c.dp;
                tmp_ado[g] = c.ad_other < 0 ? bcf_int32_missing : c.ad_other;
                // ref allele depth goes first
                tmp_ad[(var.variation.size() + 1)*g] = c.ad_ref >= 0 ? c.ad_ref : bcf_int32_missing;
                bool any_missing_gt = tmp_ad[(var.variation.size() + 1)*g] == bcf_int32_missing || var.calls[g].ngt == 0;

                // fill in other ADs if known
                for (size_t j = 0; j < var.calls[g].ngt; ++j)
                {
                    int this_ad = c.ad[j] >= 0 ? c.ad[j] : bcf_int32_missing;
                    if(   var.calls[g].gt[j] > 0
                      || (var.calls[g].gt[j] == 0 && tmp_ad[(var.variation.size() + 1)*g] == bcf_int32_missing))
                    {
                        tmp_ad[var.calls[g].gt[j] + (var.variation.size() + 1)*g] = this_ad;
                    }
                    else
                    {
                        if(var.calls[g].gt[j] < 0)
                        {
                            any_missing_gt = true;
                            break;
                        }
                    }
                }
                if(any_missing_gt)
                {
                    for(size_t qq = 0; qq < var.variation.size() + 1; ++qq)
                    {
                        tmp_ad[(var.variation.size() + 1)*g + qq] = bcf_int32_missing;
                    }
                }

                /** transfer all other formats */
                if(var.calls[g].bcf_hdr && var.calls[g].bcf_rec)
                {
                    bcf_unpack(var.calls[g].bcf_rec.get(), BCF_UN_FMT);
                    const std::set<int> skip = {
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "GT"),
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "DP"),
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "AD"),
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "ADO"),
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "AGT"),
                        // don't translate these from source bcf, they might have changed
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "AN"),
                        bcf_hdr_id2int(var.calls[g].bcf_hdr.get(), BCF_DT_ID, "AC"),
                    };
                    for(int f = 0;  f < var.calls[g].bcf_rec->n_fmt; ++f)
                    {
                        const bcf_fmt_t * fmt = &(var.calls[g].bcf_rec->d.fmt[f]);
                        if(skip.count(fmt->id))
                        {
                            continue;
                        }
                        const char * id = bcf_hdr_int2id(var.calls[g].bcf_hdr.get(), BCF_DT_ID, fmt->id);
                        switch(fmt->type)
                        {
                            case BCF_BT_INT8:
                            case BCF_BT_INT16:
                            case BCF_BT_INT32:
                            {
                                int defaultint = bcf_int32_missing;
                                int value = bcfhelpers::getFormatInt(var.calls[g].bcf_hdr.get(),
                                                                     var.calls[g].bcf_rec.get(),
                                                                     id,
                                                                     var.calls[g].bcf_sample,
                                                                     defaultint);
                                auto pf = int_fmts.find(id);
                                if(pf == int_fmts.end())
                                {
                                    pf = int_fmts.emplace(std::make_pair(
                                            std::string(id),
                                            std::vector<int>(var.calls.size(), defaultint))
                                    ).first;
                                }
                                pf->second[g] = value;
                                break;
                            }
                            case BCF_BT_FLOAT:
                            {
                                float defaultfloat;
                                memcpy(&defaultfloat, &bcf_float_missing, sizeof(float));
                                float value = (float) bcfhelpers::getFormatDouble(var.calls[g].bcf_hdr.get(),
                                                                                  var.calls[g].bcf_rec.get(),
                                                                                  id,
                                                                                  var.calls[g].bcf_sample);
                                if(std::isnan(value))
                                {
                                    memcpy(&value, &bcf_float_missing, sizeof(float));
                                }
                                auto pf = float_fmts.find(id);
                                if(pf == float_fmts.end())
                                {
                                    pf = float_fmts.emplace(std::make_pair(
                                            std::string(id),
                                            std::vector<float>(var.calls.size(), defaultfloat))
                                    ).first;
                                }
                                pf->second[g] = value;
                                break;
                            }
                                break;
                            case BCF_BT_CHAR:
                            {
                                std::string value = bcfhelpers::getFormatString(var.calls[g].bcf_hdr.get(),
                                                                                var.calls[g].bcf_rec.get(),
                                                                                id,
                                                                                var.calls[g].bcf_sample);
                                auto pf = string_fmts.find(id);
                                if(pf == string_fmts.end())
                                {
                                    pf = string_fmts.emplace(std::make_pair(
                                            std::string(id),
                                            std::vector<std::string>(var.calls.size()))
                                    ).first;
                                }
                                pf->second[g] = value;
                                break;
                            }
                            case BCF_BT_NULL:
                            default:
                                break;
                        }
                    }
                }
            }

            bcf_update_format_int32(hdr, rec, "AD", tmp_ad.get(), bcf_hdr_nsamples(hdr)*((int)var.variation.size() + 1));
            bcf_update_format_int32(hdr, rec, "ADO", tmp_ado.get(), bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "DP", tmp_dp.get(), bcf_hdr_nsamples(hdr));

            for(auto & f : int_fmts)
            {
                bcfhelpers::setFormatInts(hdr, rec, f.first.c_str(), f.second);
            }
            for(auto & f : float_fmts)
            {
                bcfhelpers::setFormatFloats(hdr, rec, f.first.c_str(), f.second);
            }
            for(auto & f : string_fmts)
            {
                bcfhelpers::setFormatStrings(hdr, rec, f.first.c_str(), f.second);
            }
        }

        // write ambiguous
        bool any_ambiguous = false;
        for (auto const & x : var.ambiguous_alleles)
        {
            if(!x.empty())
            {
                any_ambiguous = true;
            }
        }

        if(any_ambiguous)
        {
            int nsamples = bcf_hdr_nsamples(hdr);
            char ** ambiguous = new char * [nsamples];
            for (int i = 0; i < nsamples; ++i)
            {
                if (((size_t)i) >= var.ambiguous_alleles.size() || var.ambiguous_alleles[i].empty())
                {
                    ambiguous[i] = new char[2];
                    ambiguous[i][0] = '.';
                    ambiguous[i][1] = 0;
                }
                else
                {
                    std::ostringstream oss;
                    for (auto it = var.ambiguous_alleles[i].begin(); it != var.ambiguous_alleles[i].end(); ++it)
                    {
                        oss << *it;
                        if (std::next(it) != var.ambiguous_alleles[i].end())
                        {
                            oss << ",";
                        }
                    }
                    size_t N = oss.str().size()+1;
                    ambiguous[i] = new char[N];
                    strncpy(ambiguous[i], oss.str().c_str(), N);
                }
            }
            bcf_update_format_string(hdr, rec, "AGT", (const char**)ambiguous, nsamples);
            for (int i = 0; i < nsamples; ++i)
            {
                delete [] ambiguous[i];
            }
            delete [] ambiguous;
        }

        bcf_write1(_impl->fp, hdr, rec);
    }

} // namespace variant
