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
 * \brief Count variant numbers from various intermediate data structures
 *
 * \file VariantStatistics.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Variant.hh"
#include "Fasta.hh"
#include "Alignment.hh"

#include "Error.hh"

#include <memory>
#include <bitset>

// #define DEBUG_VARIANTSTATISTICS

namespace variant
{

/* variant types in uint64 */

/* encode allele types seen in low 4 bits, highest bit is unused */
static const uint64_t VT_NOCALL         = 0;
static const uint64_t VT_SNP            = 1;
static const uint64_t VT_INS            = 2;
static const uint64_t VT_DEL            = 4;
static const uint64_t VT_REF            = 8;

static const char * VT_NAMES [] = {
    "nc",               //    0x00
    "s",                //    0x01
    "i",                //    0x02
    "si",               //    0x03
    "d",                //    0x04
    "sd",               //    0x05
    "id",               //    0x06
    "sid",              //    0x07
    "r",                //    0x08
    "rs",               //    0x09
    "ri",               //    0x0a
    "rsi",              //    0x0b
    "rd",               //    0x0c
    "rsd",              //    0x0d
    "rid",              //    0x0e
    "rsid",             //    0x0f
};

/** encode call type / zygosity in bits 5-8 */
static const uint64_t CT_NUCLEOTIDES = 0;     // count nucleotides
static const uint64_t CT_ALLELES = 0x10;      // count alleles
static const uint64_t CT_HOMREF = 0x20;       // locations with only one allele seen with copy number 2
static const uint64_t CT_HET = 0x30;          // locations with one ref and one alt allele seen with copy number 1
static const uint64_t CT_HETALT = 0x40;       // locations with two alt alleles seen with copy number 1 each
static const uint64_t CT_HEMI = 0x50;         // locations with only one allele seen with copy number 1
static const uint64_t CT_AMBI = 0x60;         // locations with more than two alleles
static const uint64_t CT_HALFCALL = 0x70;     // locations with one allele, and one missing call
static const uint64_t CT_NOCALL = 0x80;       // locations with one allele, and one missing call
static const uint64_t CT_HOMALT = 0x90;       // locations with one allele, and one missing call
static const uint64_t CT_UNKNOWN = 0xa0;      // unknown locations / counts

static const char * CT_NAMES [] = {
    "nuc", // 0x00
    "al",     // 0x10
    "homref",         // 0x20
    "het",         // 0x30
    "hetalt",      // 0x40
    "hemi",        // 0x50
    "ambi",        // 0x60
    "halfcall",    // 0x70
    "nocall",      // 0x80
    "homalt",      // 0x90
    "unknown",     // 0xa0
    "unknown",     // 0xb0
    "unknown",     // 0xc0
    "unknown",     // 0xd0
    "unknown",     // 0xe0
    "unknown",     // 0xf0
};

/** we count 256 distinct things */
static const int VS_COUNTS = 256;

struct VariantStatisticsImpl
{
    VariantStatisticsImpl(const char * ref_fasta, bool _count_homref) : ref(ref_fasta), count_homref(_count_homref),
        alignment(makeAlignment("klibg"))
    {
        memset(counts, 0, sizeof(size_t)*VS_COUNTS);
        memset(rtypes, 0, sizeof(int)*VS_COUNTS);
        nrtypes = 0;
        rtype_bs.reset();
    }

    ~VariantStatisticsImpl() {}

    VariantStatisticsImpl(VariantStatisticsImpl const & rhs) : ref(rhs.ref),
        alignment(makeAlignment("klibg"))
    {
        memcpy(counts, rhs.counts, sizeof(size_t)*VS_COUNTS);
        count_homref = rhs.count_homref;
        nrtypes = rhs.nrtypes;
        memcpy(rtypes, rhs.rtypes, sizeof(int)*256);
        rtype_bs = rhs.rtype_bs;
    }

    /** translate count key to name */
    static inline std::string c2n(uint64_t _c) {
        uint8_t c = (uint8_t)_c;
        return std::string(CT_NAMES[c >> 4]) + "__" +
               std::string(VT_NAMES[c & 0x0f]);
    }

    /** translate name to count key */
    static inline uint8_t n2c(const char * n) {
        static struct _n2clookup {
            _n2clookup() {
                for(int j = 0; j < VS_COUNTS; ++j) {
                    t[c2n(j)] = j;
                }
            }
            std::map<std::string, uint8_t> t;
        } lookup;
        auto it = lookup.t.find(std::string(n));
        if(it != lookup.t.end()) {
            return it->second;
        } else {
            error("Unknown statistics key: %s", n);
        }
        return 0;
    }

    /** add single allele */
    int add_al(const char * chr, RefVar const & rhs)
    {
        size_t total_snp = 0;
        size_t total_del = 0;
        size_t total_ins = 0;
        size_t total_hom = 0;
        realignRefVar(ref, chr, rhs, alignment.get(),
                total_snp,
                total_ins,
                total_del,
                total_hom);

        // count nucleotides
        count( CT_NUCLEOTIDES | VT_SNP, total_snp );
        count( CT_NUCLEOTIDES | VT_INS, total_ins );
        count( CT_NUCLEOTIDES | VT_DEL, total_del );
        if(count_homref) {
            count( CT_NUCLEOTIDES | VT_REF, total_hom );
        }

        uint8_t t = 0;
        if(total_snp) t |= VT_SNP;
        if(total_ins) t |= VT_INS;
        if(total_del) t |= VT_DEL;

        count(CT_ALLELES | t, 1);
        return t;
    }

    size_t counts[VS_COUNTS];

    // keep track of returned types
    void count(int rt, size_t n) {
        if(!n) {
            return;
        }
        if(!rtype_bs[rt & 0xff]) {
            rtype_bs[rt & 0xff] = 1;
            rtypes[nrtypes++] = rt;
        }
        counts[rt & 0xff] += n;
    }

    void reset_rtypes() {
        nrtypes = 0;
        rtype_bs.reset();
    }

    int rtypes[256];
    int nrtypes;
    std::bitset<256> rtype_bs;

    FastaFile ref;
    bool count_homref;
    std::unique_ptr<Alignment> alignment;
};

VariantStatistics::VariantStatistics(const char * ref_fasta, bool count_homref)
{
    _impl = new VariantStatisticsImpl(ref_fasta, count_homref);
}

VariantStatistics::VariantStatistics(VariantStatistics const & rhs)
{
    _impl = new VariantStatisticsImpl(*rhs._impl);
}

VariantStatistics::~VariantStatistics()
{
    delete _impl;
}

VariantStatistics const & VariantStatistics::operator=(VariantStatistics const & rhs)
{
    delete _impl;
    _impl = new VariantStatisticsImpl(*rhs._impl);
    return *this;
}

/**
 * @brief Read and write from/to JSON files
 */
void VariantStatistics::read(Json::Value const & root)
{
    for (int i = 0; i < VS_COUNTS; ++i)
    {
        std::string k = VariantStatisticsImpl::c2n(i);
        if (root.isMember(k))
        {
            _impl->counts[i] = root[k].asUInt64();
        }
        else
        {
            _impl->counts[i] = 0;
        }
    }
}

Json::Value VariantStatistics::write()
{
    Json::Value root;
    for (int i = 0; i < VS_COUNTS; ++i)
    {
        if (_impl->counts[i])
        {
            root[VariantStatisticsImpl::c2n(i)] = Json::Value::UInt64(_impl->counts[i]);
        }
    }
#ifdef DEBUG_VARIANTSTATISTICS
    Json::StyledWriter sw;
    std::cerr << sw.write(root) << std::endl;
#endif
    return root;
}

/**
 * @brief add counts
 */
void VariantStatistics::add(VariantStatistics const & rhs)
{
    for (int i = 0; i < VS_COUNTS; ++i)
    {
        _impl->counts[i] += rhs._impl->counts[i];
    }
}

void VariantStatistics::add(Variants const & rhs, int sample, int ** rtypes, int * nrtypes)
{
    if(rtypes && nrtypes) { _impl->reset_rtypes(); }

    int types = 0;
    // count hom-alt alleles only once
    if(rhs.calls[sample].isHomalt())
    {
        types = _impl->add_al(rhs.chr.c_str(), rhs.variation[rhs.calls[sample].gt[0]-1]);
    } else if(rhs.calls[sample].isHomref()) {
        types = VT_REF;
    } else {
        for(size_t i = 0; i < rhs.calls[sample].ngt; ++i)
        {
            if (rhs.calls[sample].gt[i] > 0)
            {
                types |= _impl->add_al(rhs.chr.c_str(), rhs.variation[rhs.calls[sample].gt[i]-1]);
            }
            else if(rhs.calls[sample].gt[i] == 0)
            {
                types |= VT_REF;
            }
        }
    }

    bool is_ambi = false;
    for(auto i : rhs.ambiguous_alleles[sample])
    {
        is_ambi = true;
        if (i > 0)
        {
            types |= _impl->add_al(rhs.chr.c_str(), rhs.variation[i]);
        } else if(i == 0) {
            types |= VT_REF;
        } else {
            types |= VT_NOCALL;
        }
    }

    int location_type = CT_UNKNOWN;
    if(is_ambi) {
        location_type = CT_AMBI;
    } else if(rhs.calls[sample].isHomref()) {
        location_type = CT_HOMREF;
    } else if(rhs.calls[sample].isHet()) {
        location_type = CT_HET;
    } else if(rhs.calls[sample].isHomalt()) {
        location_type = CT_HOMALT;
    } else if(rhs.calls[sample].isHetAlt()) {
        location_type = CT_HETALT;
    } else if(rhs.calls[sample].isHalfcall()) {
        location_type = CT_HALFCALL;
    } else if(rhs.calls[sample].isNocall()) {
        location_type = CT_NOCALL;
    } else if(rhs.calls[sample].isHemi()) {
        location_type = CT_HEMI;
    }

    _impl->count(location_type | types, 1);

    if(rtypes && nrtypes) { *rtypes = _impl->rtypes; *nrtypes = _impl->nrtypes; }
}

void VariantStatistics::add(const char * chr, RefVar const & rhs, int ** rtypes, int * nrtypes)
{
    if(rtypes && nrtypes) { _impl->reset_rtypes(); }
    _impl->add_al(chr, rhs);
    if(rtypes && nrtypes) { *rtypes = _impl->rtypes; *nrtypes = _impl->nrtypes; }
}

/** resolve types to strings */
std::string VariantStatistics::type2string(int type)
{
    return VariantStatisticsImpl::c2n(type);
}

/** resolve types to strings */
int VariantStatistics::string2type(const char * str)
{
    return VariantStatisticsImpl::n2c(str);
}

}
