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

    VariantWriter & VariantWriter::operator=(VariantWriter const & rhs)
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
            stringutil::split(std::string(hdr_text, (unsigned long) len), result, "\n");
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
                // no ID => pass on
                bcf_hdr_add_hrec(hdr, hrec);
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
                    idx = bcf_hrec_find_key(hrec, "Description");
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
        rec->pos = (int32_t) var.pos;

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
            const float cqual = c.qual;
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
            for (int j = 0; j < std::max(MAX_GT, (int ) var.calls[g].ngt); ++j)
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
            // END has been written above already
            std::set<std::string> dont_write = {"AN", "AC", "END"};
            for(auto const & id : var.infos.getMemberNames())
            {
                if(dont_write.count(id))
                {
                    continue;
                }
                Json::Value const & v = var.infos[id];

                if(v.isArray() && !v.empty())
                {
                    if(v[0].type() == Json::ValueType::intValue)
                    {
                        std::unique_ptr<int[]> p_values(new int[v.size()]);
                        for(int s = 0; s < (int)v.size(); ++s)
                        {
                            try
                            {
                                p_values.get()[s] = v[s].asInt();
                            }
                            catch(std::runtime_error const & )
                            {
                                p_values.get()[s] = bcf_int32_missing;
                            }
                        }
                        bcf_update_info_int32(hdr, rec, id.c_str(), p_values.get(), (int)v.size());
                    }
                    else if(v[0].isNumeric())
                    {
                        std::unique_ptr<float[]> p_values(new float[v.size()]);
                        for(int s = 0; s < (int)v.size(); ++s)
                        {
                            try
                            {
                                p_values.get()[s] = v[s].asFloat();
                            }
                            catch(std::runtime_error const & )
                            {
                                p_values.get()[s] = bcfhelpers::missing_float();
                            }
                        }
                        bcf_update_info_float(hdr, rec, id.c_str(), p_values.get(), (int)v.size());
                    }
                    else
                    {
                        std::string value;
                        for(int s = 0; s < (int)v.size(); ++s)
                        {
                            if(s > 0)
                            {
                                value += ",";
                            }
                            value += v[s].asString();
                        }
                        bcf_update_info_string(hdr, rec, id.c_str(), value.c_str());
                    }
                }
                else if(v.isInt())
                {
                    int value;
                    try
                    {
                        value = v.asInt();
                    }
                    catch(std::runtime_error const & )
                    {
                        value = bcf_int32_missing;
                    }
                    bcf_update_info_int32(hdr, rec, id.c_str(), &value, 1);
                }
                else if(v.isNumeric())
                {
                    float value;
                    try
                    {
                        value = v.asFloat();
                    }
                    catch(std::runtime_error const & )
                    {
                        value = bcfhelpers::missing_float();
                    }
                    bcf_update_info_float(hdr, rec, id.c_str(), &value, 1);
                }
                else if(v.isBool())
                {
                    if(v.asBool())
                    {
                        bcf_update_info_flag(hdr, rec, id.c_str(), NULL, 1);
                    }
                    else
                    {
                        bcf_update_info_flag(hdr, rec, id.c_str(), NULL, 0);
                    }
                }
                else
                {
                    bcf_update_info_string(hdr, rec, id.c_str(), v.asCString());
                }
            }
        }

        if (_impl->write_formats)
        {
            std::map<std::string, std::vector<int>> int_fmts;
            std::map<std::string, std::vector<float>> float_fmts;
            std::map<std::string, std::vector<std::string>> string_fmts;

            auto tmp_dp = std::unique_ptr<int32_t[]>(new int32_t [bcf_hdr_nsamples(hdr)]);
            auto tmp_ad = std::unique_ptr<int32_t[]>(new int32_t [bcf_hdr_nsamples(hdr) * (var.variation.size() + 1)]);
            auto tmp_ado = std::unique_ptr<int32_t[]>(new int32_t [bcf_hdr_nsamples(hdr)]);

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
                        if(var.calls[g].gt[j] < (int)var.variation.size() + 1)
                        {
                            tmp_ad[var.calls[g].gt[j] + (var.variation.size() + 1)*g] = this_ad;
                        }
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
                for(auto const & id : c.formats.getMemberNames())
                {
                    Json::Value const & v = c.formats[id];
                    if(v.isArray() && !v.empty())
                    {
                        if(v[0].isInt())
                        {
                            auto f_it = int_fmts.find(id);
                            int v_size = v.size();
                            if(f_it == int_fmts.end())
                            {
                                f_it = int_fmts.emplace(std::make_pair(std::string(id),
                                                                       std::vector<int>(v.size()*var.calls.size()))).first;
                                std::fill(f_it->second.begin(), f_it->second.end(), bcf_int32_missing);
                            }
                            else if(v.size()*var.calls.size() != f_it->second.size())
                            {
                                std::cerr << "[W] Inconsistent format field counts at " << var.chr << ":" << var.pos << "\n";
                                v_size = std::min(v_size, (int) (f_it->second.size() / var.calls.size()));
                            }
                            for(auto s = 0; s < v_size; ++s)
                            {
                                try
                                {
                                    f_it->second[v_size*g + s] = v[s].asInt();
                                }
                                catch(std::runtime_error const & )
                                {
                                    f_it->second[v_size*g + s] = bcf_int32_missing;
                                }
                            }
                        }
                        else if(v[0].isNumeric())
                        {
                            auto f_it = float_fmts.find(id);
                            int v_size = v.size();
                            if(f_it == float_fmts.end())
                            {
                                f_it = float_fmts.emplace(std::make_pair(std::string(id),
                                                                         std::vector<float>(v.size()*var.calls.size()))).first;
                                union { float f; uint32_t i; } missing;
                                missing.i = bcf_float_missing;
                                std::fill(f_it->second.begin(), f_it->second.end(), missing.f);
                            }
                            else if(v.size()*var.calls.size() != f_it->second.size())
                            {
                                std::cerr << "[W] Inconsistent format field counts at " << var.chr << ":" << var.pos << "\n";
                                v_size = std::min(v_size, (int) (f_it->second.size() / var.calls.size()));
                            }
                            for(auto s = 0; s < v_size; ++s)
                            {
                                f_it->second[v_size*g + s] = v[s].asFloat();
                            }
                        }
                        else
                        {
                            auto f_it = string_fmts.find(id);
                            if(f_it == string_fmts.end())
                            {
                                f_it = string_fmts.emplace(std::make_pair(std::string(id),
                                                                          std::vector<std::string>(var.calls.size()))).first;
                            }
                            f_it->second[g] = v.asString();
                        }
                    }
                    else if(v.isInt())
                    {
                        auto f_it = int_fmts.find(id);
                        int v_size = 1;
                        if(f_it == int_fmts.end())
                        {
                            f_it = int_fmts.emplace(std::make_pair(std::string(id),
                                                                   std::vector<int>(var.calls.size()))).first;
                            std::fill(f_it->second.begin(), f_it->second.end(), bcf_int32_missing);
                        }
                        else if(var.calls.size() != f_it->second.size())
                        {
                            std::cerr << "[W] Inconsistent format field counts at " << var.chr << ":" << var.pos << "\n";
                            v_size = (int) (f_it->second.size() / var.calls.size());
                        }
                        try
                        {
                            f_it->second[v_size*g] = v.asInt();
                        }
                        catch(std::runtime_error const & )
                        {
                            f_it->second[v_size*g] = bcf_int32_missing;
                        }
                    }
                    else if(v.isNumeric())
                    {
                        auto f_it = float_fmts.find(id);
                        int v_size = 1;
                        if(f_it == float_fmts.end())
                        {
                            f_it = float_fmts.emplace(std::make_pair(std::string(id),
                                                                     std::vector<float>(var.calls.size()))).first;
                            union { float f; uint32_t i; } missing;
                            missing.i = bcf_float_missing;
                            std::fill(f_it->second.begin(), f_it->second.end(), missing.f);
                        }
                        else if(var.calls.size() != f_it->second.size())
                        {
                            std::cerr << "[W] Inconsistent format field counts at " << var.chr << ":" << var.pos << "\n";
                            v_size = (int) (f_it->second.size() / var.calls.size());
                        }
                        f_it->second[v_size*g] = v.asFloat();
                    }
                    else
                    {
                        auto f_it = string_fmts.find(id);
                        if(f_it == string_fmts.end())
                        {
                            f_it = string_fmts.emplace(std::make_pair(std::string(id),
                                                                      std::vector<std::string>(var.calls.size()))).first;
                        }
                        f_it->second[g] = v.asString();
                    }
                }
            }

            bcf_update_format_int32(hdr, rec, "AD", tmp_ad.get(), bcf_hdr_nsamples(hdr)*((int)var.variation.size() + 1));
            bcf_update_format_int32(hdr, rec, "ADO", tmp_ado.get(), bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "DP", tmp_dp.get(), bcf_hdr_nsamples(hdr));

            for(auto & i : int_fmts)
            {
                bcf_update_format_int32(hdr, rec, i.first.c_str(), i.second.data(), (int)i.second.size());
            }
            for(auto & f : float_fmts)
            {
                bcf_update_format_float(hdr, rec, f.first.c_str(), f.second.data(), (int)f.second.size());
            }
            for(auto & s : string_fmts)
            {
                std::unique_ptr<const char *[]> values(new const char *[s.second.size()]);
                for(int ix = 0; ix < (int)s.second.size(); ++ix)
                {
                    values.get()[ix] = s.second[ix].c_str();
                }
                bcf_update_format_string(hdr, rec, s.first.c_str(), values.get(), (int)s.second.size());
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
