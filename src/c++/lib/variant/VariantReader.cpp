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
 * \file Variant.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include "VariantImpl.hh"

//#define DEBUG_VARIANT_GTS

namespace variant
{

VariantReader::VariantReader()
{
    _impl = new VariantReaderImpl();
}

VariantReader::VariantReader(VariantReader const & rhs)
{
    _impl = new VariantReaderImpl();
    for (size_t i = 0; i < rhs._impl->samples.size(); ++i)
    {
        addSample(rhs._impl->samples[i].filename.c_str(),
                  rhs._impl->samples[i].sample.c_str());
    }
    if (rhs._impl->regions != "")
    {
        setRegions(rhs._impl->regions.c_str(), rhs._impl->regionsFile);
    }
    if (rhs._impl->targets != "")
    {
        setTargets(rhs._impl->targets.c_str(), rhs._impl->targetsFile);
    }
    _impl->applyFilters = rhs._impl->applyFilters;
    _impl->applyFiltersPerSample = rhs._impl->applyFiltersPerSample;
    _impl->returnHomref = rhs._impl->returnHomref;
    _impl->buffered_variants = rhs._impl->buffered_variants;
}

VariantReader::~VariantReader()
{
    delete _impl;
}

VariantReader const & VariantReader::operator=(VariantReader const & rhs)
{
    if(&rhs == this)
    {
        return *this;
    }
    delete _impl;

    _impl = new VariantReaderImpl();
    for (size_t i = 0; i < rhs._impl->samples.size(); ++i)
    {
        addSample(rhs._impl->samples[i].filename.c_str(),
                  rhs._impl->samples[i].sample.c_str());
    }
    if (rhs._impl->regions != "")
    {
        setRegions(rhs._impl->regions.c_str(), rhs._impl->regionsFile);
    }
    if (rhs._impl->targets != "")
    {
        setTargets(rhs._impl->targets.c_str(), rhs._impl->targetsFile);
    }
    _impl->applyFilters = rhs._impl->applyFilters;
    _impl->applyFiltersPerSample = rhs._impl->applyFiltersPerSample;
    _impl->returnHomref = rhs._impl->returnHomref;
    _impl->buffered_variants = rhs._impl->buffered_variants;
    return *this;
}

void VariantReader::setApplyFilters(bool filters, int sample)
{
    if(sample < 0)
    {
        _impl->applyFilters = filters;
    }
    else
    {
        if(int(_impl->applyFiltersPerSample.size()) <= sample)
        {
            _impl->applyFiltersPerSample.resize((unsigned long) (sample + 1), false);
            _impl->applyFiltersPerSample[sample] = filters;
        }
    }
}

/**
 * @brief change GTs on chrX/Y to be diploid for matching
 *
 */
void VariantReader::setFixChrXGTs(bool fix)
{
    _impl->fix_chrX = fix;
}

bool VariantReader::getApplyFilters(int sample) const
{
    if(sample < 0)
    {
        return _impl->applyFilters;
    }
    else
    {
        return _impl->applyFilters || (int(_impl->applyFiltersPerSample.size()) > sample && _impl->applyFiltersPerSample[sample]);
    }
}

/**
 * @brief Return homref/no-calls
 *
 */
void VariantReader::setReturnHomref(bool homref)
{
    _impl->returnHomref = homref;
}

bool VariantReader::getReturnHomref() const
{
    return _impl->returnHomref;
}

/**
 * @brief Validate reference alleles
 *
 */
void VariantReader::setValidateRef(const char * ref_fasta, bool validate)
{
    _impl->ref = FastaFile(ref_fasta);
    _impl->validateRef = validate;
}

bool VariantReader::getValidateRef() const
{
    return _impl->validateRef;
}

/**
 * @brief Add a sample to read from
 * @param filename  file name
 * @param sname  name of sample to read from
 * @return sample number to retrieve records in get()
 */
int VariantReader::addSample(const char * filename, const char * sname)
{
    SampleInfo si;

    si.filename = filename;
    si.sample = sname;

    si.ireader = -1;
    auto psi = _impl->filename_mapping.find(filename);
    if(psi != _impl->filename_mapping.end())
    {
        si.ireader = psi->second;
    }
    else
    {
        if (!bcf_sr_add_reader(_impl->files, filename))
        {
            error("Failed to open or file not indexed: %s\n", filename);
        }
        si.ireader = _impl->files->nreaders - 1;
        _impl->filename_mapping[filename] = si.ireader;
        bcf_hdr_append(_impl->files->readers->header,
                       "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that could not be imported.\">");
        bcf_hdr_sync(_impl->files->readers->header);
    }

    if(sname && strlen(sname) > 0)
    {
        if (std::string(sname) == "*")
        {
            int rv = (int)_impl->samples.size();

            for (int i = 0; i < bcf_hdr_nsamples(_impl->files->readers[si.ireader].header); ++i)
            {
                if(std::string(_impl->files->readers[si.ireader].header->samples[i]) == "*")
                {
                    std::cerr << "Skipping sample named '*'" << "\n";
                    continue;
                }
                /* std::cerr << "Adding sample: " << _impl->files->readers[si.ireader].header->samples[i] << "\n"; */
                addSample(filename, _impl->files->readers[si.ireader].header->samples[i]);
            }

            return rv;
        }
        else
        {
            si.isample = bcf_hdr_id2int(_impl->files->readers[si.ireader].header,
                                        BCF_DT_SAMPLE,
                                        sname);

            if(si.isample < 0)
            {
                error("Sample %s doesn't exist in file %s\n", sname, filename);
            }
            si.end_in_vcf = -1;
            _impl->samples.push_back(si);
        }
    }
    else
    {
        si.isample = 0;
        si.end_in_vcf = -1;
        _impl->samples.push_back(si);
    }
    return (int)(_impl->samples.size() - 1);
}

void VariantReader::getSampleList(std::list< std::pair<std::string, std::string> > & files)
{
    for (SampleInfo const & si : _impl->samples)
    {
        files.push_back(std::make_pair(si.filename, si.sample));
    }
}


/**
 * @brief Interface to htslib regions functionality
 * @param regions regions string, see synced_bcf_reader.h
 * @param isFile True if regions is a file
 */
void VariantReader::setRegions(const char * regions, bool isFile)
{
    int result = bcf_sr_set_regions(_impl->files,
                                    regions, isFile ? 1 : 0);
    if(result < 0)
    {
        error("Failed to set regions string %s.", regions);
    }
    _impl->regions = regions;
    _impl->regionsFile = isFile;
}

/**
 * @brief Interface to htslib targets functionality
 * @param targets targets string, see synced_bcf_reader.h
 * @param isFile True if targets is a file
 */
void VariantReader::setTargets(const char * targets, bool isFile)
{
    int result = bcf_sr_set_targets(_impl->files,
                                    targets, isFile ? 1 : 0, 0);
    if(result < 0)
    {
        error("Failed to set targets string %s.", targets);
    }
    _impl->targets = targets;
    _impl->targetsFile = isFile;
}


/**
 * @brief Rewind / set region to read
 *
 * @param chr chromosome/contig name
 * @param startpos start position on chr (or -1 for entire chr)
 */
void VariantReader::rewind(const char * chr, int64_t startpos)
{
    int returned = -1;

    if(!chr || strcmp(chr, "") == 0)
    {
        chr = NULL;
    }

    if(startpos < 0)
    {
        returned = bcf_sr_seek(_impl->files, chr, 0);
    }
    else
    {
        returned = bcf_sr_seek(_impl->files, chr, (int) startpos);
    }
    if (returned == -_impl->files->nreaders)
    {
        error("Could not seek to: %s", stringutil::formatPos(chr, startpos).c_str());
    }
    _impl->buffered_variants.clear();
}

/**
 * @brief Return next variant and advance
 *
 * @param v Variant record to populate
 */
Variants & VariantReader::current()
{
    if (_impl->buffered_variants.empty())
    {
        static Variants tmp;
        return tmp;
    }
    else
    {
        return _impl->buffered_variants.front();
    }
}

bool VariantReader::advance()
{
    if (!_impl->buffered_variants.empty())
    {
        _impl->buffered_variants.pop_front();
    }

    if (!_impl->buffered_variants.empty())
    {
        return true;
    }

    int nl = bcf_sr_next_line(_impl->files);
    if (nl <= 0)
    {
        return false;
    }

    Variants vars;

    std::map<std::string, int> vl;

    vars.chr = "";
    vars.variation.clear();
    // calculate pos and len
    vars.pos = 0;
    vars.len = 0;
    vars.ambiguous_alleles.clear();
    bool import_fail = false;
    for (auto & si : _impl->samples)
    {
        bcf_sr_t & reader(_impl->files->readers[si.ireader]);
        if(!bcf_sr_has_line(_impl->files, si.ireader))
        {
            continue;
        }
        bcf1_t *line = reader.buffer[0];
        bcf_unpack(line, BCF_UN_SHR);

        if(vars.chr == "")
        {
            vars.chr = reader.header->id[BCF_DT_CTG][line->rid].key;
        }
        else
        {
            std::string chr2 = reader.header->id[BCF_DT_CTG][line->rid].key;
            if (vars.chr != chr2)
            {
                error("Chromosome mismatch: %s != %s", vars.chr.c_str(), chr2.c_str());
            }
        }

        int64_t refstart = line->pos;
        int64_t refend = refstart;

        int endfield = bcfhelpers::getInfoInt(reader.header, line, "END", -1);

        if(endfield > 0)
        {
            // if there is an end field, don't validate the ref allele
            refend = endfield-1;
        }
        else
        {
            refend = refstart + strlen(line->d.allele[0]) - 1;
            if (strchr(line->d.allele[0], '.') || strchr(line->d.allele[0], '-') || strchr(line->d.allele[0], '*'))
            {
                // length might be inaccurate now
                refend = refstart;
                std::cerr << "[W] Unsupported REF allele with undefined length at " << vars.chr << ":" << refstart << "-" << refend
                          << " REF allele: " << line->d.allele[0] << "\n";
                if(!import_fail)
                {
                    vars.setInfo("IMPORT_FAIL", true);
                    import_fail = true;
                }
            }
            if(_impl->validateRef)
            {
                std::string actual_ref_seq = _impl->ref.query(vars.chr.c_str(), refstart, refend);
                if(actual_ref_seq != line->d.allele[0]) {
                    std::cerr << "[W] REF allele mismatch with reference sequence at " << vars.chr << ":" << refstart << "-" << refend
                              << " VCF: " << line->d.allele[0] << " REF: " << actual_ref_seq << "\n";

                    if(!import_fail)
                    {
                        vars.setInfo("IMPORT_FAIL", true);
                        import_fail = true;
                    }
                }
            }
        }

        // use the position and length of the rightmost variant as a start
        if(vars.pos < refstart)
        {
            vars.pos = refstart;
            vars.len = refend - refstart + 1;
        }
        else if (vars.pos == refstart)
        {
            // use the maximum length over all variants starting at the final vars.pos
            vars.len = std::max(vars.len, refend - refstart + 1);
        }

        si.end_in_vcf = (int) refend;
        si.allele_map.resize(line->n_allele-1);
        for (int j = 1; j < line->n_allele; ++j)
        {
            RefVar rv;
            rv.start = refstart;
            rv.end = refend;
            rv.alt = line->d.allele[j];
            boost::to_upper(rv.alt);

            if (rv.alt == "<DEL>")
            {
                rv.alt = "";
            }
            else if(rv.alt == "*" || rv.alt == "<NON_REF>" || import_fail)
            {
                // these are OK to let through.
                // we might want to add some better matching logic later
                rv.alt = "*";
            }
            else
            {
                // check alleles
                for (size_t qq = 0; qq < rv.alt.size(); ++qq)
                {
                    /* A	A	Adenine */
                    /* C	C	Cytosine */
                    /* G	G	Guanine */
                    /* T	T	Thymine */
                    /* U	U	Uracil */
                    /* R	A or G	puRine */
                    /* Y	C, T or U	pYrimidines */
                    /* K	G, T or U	bases which are Ketones */
                    /* M	A or C	bases with aMino groups */
                    /* S	C or G	Strong interaction */
                    /* W	A, T or U	Weak interaction */
                    /* B	not A (i.e. C, G, T or U)	B comes after A */
                    /* D	not C (i.e. A, G, T or U)	D comes after C */
                    /* H	not G (i.e., A, C, T or U)	H comes after G */
                    /* V	neither T nor U (i.e. A, C or G)	V comes after U */
                    /* N	A C G T U	Nucleic acid */
                    /* X	masked */
                    /* -/.	gap of indeterminate length - not supported */
                    if(rv.alt[qq] != 'A' && \
                   rv.alt[qq] != 'C' && \
                   rv.alt[qq] != 'G' && \
                   rv.alt[qq] != 'T' && \
                   rv.alt[qq] != 'U' && \
                   rv.alt[qq] != 'R' && \
                   rv.alt[qq] != 'Y' && \
                   rv.alt[qq] != 'K' && \
                   rv.alt[qq] != 'M' && \
                   rv.alt[qq] != 'S' && \
                   rv.alt[qq] != 'W' && \
                   rv.alt[qq] != 'B' && \
                   rv.alt[qq] != 'D' && \
                   rv.alt[qq] != 'H' && \
                   rv.alt[qq] != 'V' && \
                   rv.alt[qq] != 'N' && \
                   rv.alt[qq] != 'X')
                    {
//                    std::cerr << "Invalid allele " << rv.alt.c_str() << " at " << refstart << "\n";
                        vars.setInfo("IMPORT_FAIL", true);
                        import_fail = true;
                    }
                }
            }

            // use the position and length of the rightmost variant allele as a start
            int64_t vend = vars.pos + vars.len - 1;
            if(vars.pos > rv.start)
            {
                vars.pos = rv.start;
            }

            if (vend < rv.end)
            {
                // use the maximum length over all variants starting at the final vars.pos
                vars.len = rv.end - vars.pos;
            }

            std::string rvr(rv.repr());

            if(import_fail || rv.alt == "*")
            {
                si.allele_map[j-1] = -1;
            }
            else
            {
                auto al = vl.find(rvr);
                if(al == vl.end())
                {
                    vars.variation.push_back(rv);
                    si.allele_map[j-1] = (int) vars.variation.size();
                    vl[rvr] = (int) vars.variation.size();
                }
                else
                {
                    si.allele_map[j-1] = al->second;
                }
            }
        }
    }

#ifdef DEBUG_VARIANT_GTS
    std::cerr << vars.chr << ":" << vars.pos << "\n";
#endif

    vars.calls.resize(_impl->samples.size());
    vars.ambiguous_alleles.resize(_impl->samples.size());
    for (auto & li : vars.ambiguous_alleles)
    {
        li.clear();
    }

    int ncalls = 0;
    int n_non_ref_calls = 0;

    for (size_t sid = 0; sid < _impl->samples.size(); ++sid)
    {
        SampleInfo & si = _impl->samples[sid];

        if(!bcf_sr_has_line(_impl->files, si.ireader))
        {
            vars.calls[sid].ngt = 0;
            vars.calls[sid].phased = false;
            vars.calls[sid].nfilter = 0;
            continue;
        }

        bcf_sr_t & reader(_impl->files->readers[si.ireader]);
        const int isample = si.isample;
        bcf1_t *line = reader.buffer[0];

        bcf_unpack(line, BCF_UN_FLT);

        vars.calls[sid].nfilter = (size_t) line->d.n_flt;
        if(vars.calls[sid].nfilter > MAX_FILTER)
        {
            error("Too many filters at %s:%i in sample %i", vars.chr.c_str(), vars.pos, sid);
        }

        bool fail = false;
        for(int j = 0; j < (int)vars.calls[sid].nfilter; ++j)
        {
            std::string filter = "PASS";
            int k = line->d.flt[j];

            if(k >= 0)
            {
                filter = bcf_hdr_int2id(reader.header, BCF_DT_ID, line->d.flt[j]);
            }
            if(filter != "PASS")
            {
                fail = true;
            }
            vars.calls[sid].filter[j] = filter;
        }

        if(getApplyFilters((int) sid) && fail)
        {
            vars.calls[sid].ngt = 0;
            vars.calls[sid].phased = false;
            vars.calls[sid].nfilter = 0;
            continue;
        }

        ++ncalls;
        bcf_unpack(line, BCF_UN_ALL);

        std::string fmt_strings = bcfhelpers::getFormatString(reader.header, line, "FT", isample, "");
        std::vector<std::string> fmt_filters;
        stringutil::split(fmt_strings, fmt_filters, ";", false);
        for(auto const & f : fmt_filters)
        {
            if(f.empty() || f == "PASS")
            {
                continue;
            }
            fail = true;
            bool has_filter = false;
            for(size_t ff = 0; ff < vars.calls[sid].nfilter; ++ff)
            {
                if(f == vars.calls[sid].filter[ff])
                {
                    has_filter = true;
                    break;
                }
            }
            if(!has_filter)
            {
                if(vars.calls[sid].nfilter + 1 > MAX_FILTER)
                {
                    error("Too many filters at %s:%i in sample %i", vars.chr.c_str(), vars.pos, sid);
                }
                vars.calls[sid].filter[vars.calls[sid].nfilter++] = f;
            }
        }

        if(getApplyFilters((int) sid) && fail)
        {
            vars.calls[sid].ngt = 0;
            vars.calls[sid].phased = false;
            vars.calls[sid].nfilter = 0;
            continue;
        }

        vars.calls[sid].qual = line->qual;

        for(int ni = 0; ni < line->n_info; ++ni)
        {
            bcf_info_t * inf = &line->d.info[ni];
            const char * id = bcf_hdr_int2id(reader.header, BCF_DT_ID, inf->key);
            if(vars.infos.isMember(id))
            {
                // only take first instance.
                // TODO warn or handle
                continue;
            }
            switch(inf->type)
            {
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                {
                    auto ints = bcfhelpers::getInfoInts(reader.header, line, id);
                    if(ints.size() == 1)
                    {
                        vars.infos[id] = ints[0];
                    }
                    else if(1 < ints.size())
                    {
                        vars.infos[id] = Json::Value(Json::arrayValue);
                        for(int q = 0; q < (int)ints.size(); ++q)
                        {
                            vars.infos[id][q] = ints[q];
                        }
                    }
                    break;
                }
                case BCF_BT_FLOAT:
                {
                    auto floats = bcfhelpers::getInfoFloats(reader.header, line, id);
                    if (floats.size() == 1)
                    {
                        vars.infos[id] = floats[0];
                    }
                    else if (1 < floats.size())
                    {
                        vars.infos[id] = Json::Value(Json::arrayValue);
                        for (int q = 0; q < (int)floats.size(); ++q)
                        {
                            vars.infos[id][q] = floats[q];
                        }
                    }
                    break;
                }
                case BCF_BT_CHAR:
                    vars.infos[id] = bcfhelpers::getInfoString(reader.header, line, id);
                    break;
                default:
                    vars.infos[id] = true;
                    break;
            }
        }

        int ngt = 0;
        memset(vars.calls[sid].gt, -1, MAX_GT*sizeof(int));
        bcfhelpers::getGT(reader.header, line, isample,
                          vars.calls[sid].gt,
                          ngt,
                          vars.calls[sid].phased);

        if(ngt > MAX_GT)
        {
            vars.setInfo("IMPORT_FAIL", true);
            ngt = MAX_GT;
        }

        // HAP-254: simple fix: duplicate haploid calls onto other haplotype on chrX and Y
        if(_impl->fix_chrX &&
           (vars.chr == "chrX" || vars.chr == "X" || vars.chr == "chrY" || vars.chr == "Y") &&
           ngt == 1)
        {
            ngt = 2;
            vars.calls[sid].gt[1] = vars.calls[sid].gt[0];
        }

        vars.calls[sid].ngt = (size_t) ngt;

        int adcount = std::max((int)(2*vars.calls.size()), int(vars.variation.size() + 1));
        int * ad = new int[adcount];
        memset(ad, -1, sizeof(int)*adcount);

        bcfhelpers::getAD(reader.header, line, isample,
                          ad, adcount);

        vars.calls[sid].ad_ref = ad[0];
        vars.calls[sid].ad_other = 0;
        memset(vars.calls[sid].ad, 0, sizeof(int)*MAX_GT);

        // initialize AD
        for (int j = 0; j < ngt; ++j)
        {
            if(vars.calls[sid].gt[j] >= 0 && vars.calls[sid].gt[j] < adcount
               && ad[vars.calls[sid].gt[j]] >= 0)
            {
                vars.calls[sid].ad[j] = ad[vars.calls[sid].gt[j]];
                ad[vars.calls[sid].gt[j]] = -1;
            }
        }

        for (int i = 0; i < adcount; ++i)
        {
            if(ad[i] > 0)
            {
                vars.calls[sid].ad_other += ad[i];
            }
        }

        delete [] ad;

#ifdef DEBUG_VARIANT_GTS
        std::cerr << "\ts" << sid << ": ";
        for (size_t i = 0; i < vars.calls[sid].ngt; ++i)
        {
            std::cerr << vars.calls[sid].gt[i] << " ";
        }
        std::cerr << "\n";
#endif
        // check if the GTs are valid + remap
        int gts_gte_0 = 0;
        for (size_t i = 0; i < vars.calls[sid].ngt; ++i)
        {
            int gtv = vars.calls[sid].gt[i];
            // remap all non-ref alleles
            if (gtv >= 0)
            {
                ++gts_gte_0;
            }
            if(gtv > 0)
            {
                ++n_non_ref_calls;
                if((size_t)(gtv-1) < si.allele_map.size())
                {
                    vars.calls[sid].gt[i] = si.allele_map[gtv-1];
                }
                else
                {
                    if (!import_fail)
                    {
                        vars.setInfo("IMPORT_FAIL", true);
                        import_fail = true;
                    }
                    std::cerr << "Invalid GT at " << vars.chr << ":" << vars.pos << " in sample" << sid << "\n";
                    // turn this into a no-call so it doesn't get lost later on
                    vars.calls[sid].gt[i] = -1;
                    break;
                }
            }
        }
        bcfhelpers::getDP(reader.header, line, isample,
                          vars.calls[sid].dp);

        const std::set<int> skip = {
            // these are special and have been handled above
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "GT"),
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "DP"),
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "AD"),
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "ADO"),
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "AGT"),
            // don't translate these from source bcf, they might have changed
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "AN"),
            bcf_hdr_id2int(reader.header, BCF_DT_ID, "AC"),
        };

        for(int f = 0;  f < line->n_fmt; ++f)
        {
            const bcf_fmt_t * fmt = &(line->d.fmt[f]);
            if(skip.count(fmt->id))
            {
                continue;
            }
            const char * id = bcf_hdr_int2id(reader.header, BCF_DT_ID, fmt->id);
            switch(fmt->type)
            {
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                {
                    const std::vector<int> values = bcfhelpers::getFormatInts(reader.header,
                                                                              line,
                                                                              id,
                                                                              isample);
                    if(values.empty())
                    {
                        // TODO warn?
                    }
                    if(values.size() == 1)
                    {
                        vars.calls[sid].formats[id] = values[0];
                    }
                    else
                    {
                        vars.calls[sid].formats[id] = Json::Value(Json::arrayValue);
                        for(int ffv = 0; ffv < (int)values.size(); ++ffv)
                        {
                            vars.calls[sid].formats[id] = values[ffv];
                        }
                    }
                    break;
                }
                case BCF_BT_FLOAT:
                {
                    const std::vector<float> values = bcfhelpers::getFormatFloats(reader.header,
                                                                                  line,
                                                                                  id,
                                                                                  isample);
                    if(values.empty())
                    {
                        // TODO warn?
                    }
                    if(values.size() == 1)
                    {
                        vars.calls[sid].formats[id] = values[0];
                    }
                    else
                    {
                        vars.calls[sid].formats[id] = Json::Value(Json::arrayValue);
                        for(int ffv = 0; ffv < (int)values.size(); ++ffv)
                        {
                            vars.calls[sid].formats[id] = values[ffv];
                        }
                    }
                    break;
                }
                case BCF_BT_CHAR:
                {
                    const std::string value = bcfhelpers::getFormatString(reader.header,
                                                                          line,
                                                                          id,
                                                                          isample);
                    vars.calls[sid].formats[id] = value;
                    break;
                }
                case BCF_BT_NULL:
                default:
                {
                    // TODO handle
                    break;
                }
            }
        }
    }

    // no calls unpacked because everything is filtered -> go again
    if (ncalls == 0 || ((!_impl->returnHomref) && n_non_ref_calls == 0))
    {
        return advance();
    }

    _impl->buffered_variants.push_back(vars);

    return true;
}

/**
 * @brief Insert a Variants record into the stream
 * @details The next record returned will be this one
 *
 * @param s Variants to put on top of stack
 * @param back enqueue at back or front of buffer
 */
void VariantReader::enqueue(Variants const & v, bool back)
{
    if (back)
    {
        _impl->buffered_variants.push_back(v);
    }
    else
    {
        _impl->buffered_variants.push_front(v);
    }
}


} // namespace variant

