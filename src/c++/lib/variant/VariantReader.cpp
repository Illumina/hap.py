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
#include "VariantImpl.hh"

// #define DEBUG_VARIANT_GTS

namespace variant
{

/**
 * @brief Classify a variant's GT type
 *
 */
gttype getGTType(Call const& var)
{
    if(var.ngt > 0)
    {
        if (var.ngt == 1)
        {
            if(var.gt[0] > 0)
            {
                return gt_haploid;
            }
            else if(var.gt[0] == 0)
            {
                return gt_homref;
            }
        }
        else if (var.ngt == 2)
        {
            if(var.gt[0] == 0 && var.gt[1] == 0)
            {
                return gt_homref;
            }
            else if(    (var.gt[0] == 0 && var.gt[1] > 0)
                     || (var.gt[1] == 0 && var.gt[0] > 0)  )
            {
                return gt_het;
            }
            else if(    var.gt[0] > 0 && var.gt[1] > 0
                     && var.gt[0] == var.gt[1] )
            {
                return gt_homalt;
            }
            else if(    var.gt[0] > 0 && var.gt[1] > 0
                     && var.gt[0] != var.gt[1] )
            {
                return gt_hetalt;
            }
        }
    }

    return gt_unknown;
}

std::ostream & operator<<(std::ostream &o, gttype const & v)
{
    switch(v)
    {
        case gt_haploid:
            o << "gt_haploid";
            break;
        case gt_homref:
            o << "gt_homref";
            break;
        case gt_homalt:
            o << "gt_homalt";
            break;
        case gt_het:
            o << "gt_het";
            break;
        case gt_hetalt:
            o << "gt_hetalt";
            break;
        case gt_unknown:
        default:
            o << "gt_unknown";
            break;
    }

    return o;
}

std::ostream & operator<<(std::ostream &o, Call const & v)
{
    if(v.ngt == 0)
    {
        o << ".";
    }
    for (size_t i = 0; i < v.ngt; ++i)
    {
        if(i > 0)
        {
            o << (v.phased ? "|" : "/");
        }
        o << v.gt[i];
    }

    if(v.nfilter > 0)
    {
        o << " ";
        for (size_t i = 0; i < v.nfilter; ++i)
        {
            if(i > 0)
            {
                o << ",";
            }
            o << v.filter[i];
        }
    }

    return o;
}

std::ostream & operator<<(std::ostream &o, Variants const & v)
{
    o << v.chr << ":" << v.pos << "-" << (v.pos + v.len - 1);

    for (RefVar const & rv : v.variation)
    {
        o << " " << rv;
    }

    for (Call const& c : v.calls)
    {
        o << " " << c;
    }

    bool any_ambiguous = false;
    for (auto & x : v.ambiguous_alleles)
    {
        if(!x.empty())
        {
            any_ambiguous = true;
            break;
        }
    }

    if (any_ambiguous)
    {
        o << "ambig[";
        for (auto & x : v.ambiguous_alleles)
        {
            for(auto y : x)
            {
                o << y << " ";
            }
            o << ";";
        }
        o << "]";
    }
    return o;
}

uint64_t Variants::MAX_VID = 0;
Variants::Variants() : id(MAX_VID++) {}

void setVariantInfo(Variants & v, std::string const & name, std::string const & value)
{
    size_t start = 0, p = v.info.find_first_of(name + "=");
    std::string resulting_info;
    do
    {
        if(p != std::string::npos)
        {
            size_t p2 = p+1;
            if(p == 0 || v.info[p-1] == ';')
            {
                p2 = v.info.find_first_of(";", p+1);
                resulting_info += v.info.substr(start, p-start);
                if (p2 == std::string::npos)
                {
                    start = v.info.size();
                    break;
                }
                else
                {
                    start = p2 + 1;
                }
            }
            p = v.info.find_first_of(name + "=", p2);
        }
    }
    while(p != std::string::npos);
    resulting_info += v.info.substr(start);
    if(resulting_info.size() > 0 && resulting_info[resulting_info.size()-1] == ';')
    {
        resulting_info.resize(resulting_info.size() - 1);
    }
    if(value.size() > 0)
    {
        if(resulting_info.size() > 0)
        {
            resulting_info += std::string(";");
        }
        resulting_info += name + "=" + value;
    }
    v.info = resulting_info;
}

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
    }

    si.phdr = bcfhelpers::ph(bcf_hdr_dup(_impl->files->readers[si.ireader].header));

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

bool VariantReader::advance(bool get_calls, bool get_info)
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
    vars.info = "";
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

        if(get_info)
        {
            // extract vcf info.
            // somewhat circular, when parsing a VCF, this has been read, parsed, and now we write it back
            kstring_t s;
            s.s = NULL;
            s.l = s.m = 0;
            vcf_format(reader.header, line, &s);
            size_t p0 = 0, q0 = 0, len = 0;
            int tabs = 0;
            while(tabs < 8 && p0 < s.l)
            {
                ++p0;
                if(tabs < 7)
                {
                    q0 = p0;
                }
                else
                {
                    ++len;
                }
                if (s.s[p0] == '\t')
                {
                    ++tabs;
                    if (tabs == 7)
                    {
                        ++q0;
                    }
                }
            }
            if (len > 0)
            {
                --len;
                std::string qinf(s.s + q0, len);
                if (qinf != ".")
                {
                    if (vars.info.size() > 0)
                    {
                        vars.info += ";";
                    }
                    vars.info += qinf;
                }
            }
            free(s.s);
        }

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
            if (strchr(line->d.allele[0], '.') || strchr(line->d.allele[0], '-'))
            {
                // length might be inaccurate now
                refend = refstart;
                std::cerr << "[W] Unsupported REF allele with undefined length at " << vars.chr << ":" << refstart << "-" << refend
                          << " REF allele: " << line->d.allele[0] << "\n";
                if(!import_fail)
                {
                    if (!vars.info.empty())
                    {
                        vars.info += ";";
                    }
                    vars.info += "IMPORT_FAIL";
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
                        if (!vars.info.empty())
                        {
                            vars.info += ";";
                        }
                        vars.info += "IMPORT_FAIL";
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
                continue;
            }

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
                    std::cerr << "Invalid allele " << rv.alt.c_str() << " at " << refstart << "\n";
                    if (!vars.info.empty())
                    {
                        vars.info += ";";
                    }
                    vars.info += "IMPORT_FAIL";
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

#ifdef DEBUG_VARIANT_GTS
    std::cerr << vars.chr << ":" << vars.pos << "\n";
#endif

    if(get_calls)
    {
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

            bcfhelpers::p_bcf1 p_line = bcfhelpers::pb(bcf_dup(line));

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

            int ngt = 0;
            bcfhelpers::getGT(reader.header, line, isample,
                              vars.calls[sid].gt,
                              ngt,
                              vars.calls[sid].phased);

            if(ngt > MAX_GT)
            {
                if (vars.info != "")
                {
                    vars.info += ";";
                }
                vars.info += "IMPORT_FAIL_GT_" + std::to_string(sid + 1);
                ngt = MAX_GT;
            }

            vars.calls[sid].ngt = (size_t) ngt;

            int adcount = int(vars.variation.size() + 1);
            int * ad = new int[adcount];
            memset(ad, -1, sizeof(int)*adcount);

            bcfhelpers::getAD(reader.header, line, isample,
                              ad, adcount);

            vars.calls[sid].ad_ref = ad[0];
            vars.calls[sid].ad_other = 0;

            // initialize AD
            for (int j = 0; j < ngt; ++j)
            {
                if(vars.calls[sid].gt[j] >= 0 && ad[vars.calls[sid].gt[j]] >= 0)
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
                            if (!vars.info.empty())
                            {
                                vars.info += ";";
                            }
                            vars.info += "IMPORT_FAIL";
                            import_fail = true;
                        }
                        std::cerr << "Invalid GT at " << vars.chr << ":" << vars.pos << " in sample" << sid << "\n";
                        // turn this into a homref call so it doesn't get lost later on
                        vars.calls[sid].gt[i] = 0;
                        break;
                    }
                }
            }
            // normalize no-calls
            if (gts_gte_0 == 0)
            {
                vars.calls[sid].ngt = 0;
            }
            vars.calls[sid].bcf_hdr = si.phdr;
            vars.calls[sid].bcf_rec = p_line;
            vars.calls[sid].bcf_sample = isample;

            bcfhelpers::getDP(reader.header, line, isample,
                              vars.calls[sid].dp);

        }

        // no calls unpacked because everything is filtered -> go again
        if (ncalls == 0 || ((!_impl->returnHomref) && n_non_ref_calls == 0))
        {
            return advance(get_calls, get_info);
        }
    }
    else
    {
        vars.calls.clear();
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

