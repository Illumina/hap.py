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

#include <memory>

// #define DEBUG_VARIANTSTATISTICS

namespace variant
{

/* variant types in uint64 */

/* encode allele type in 4 bits */

static const uint64_t VT_UNKNOWN        = 0;
static const uint64_t VT_SNP            = 1; // 0001b
static const uint64_t VT_INS            = 2; // 0010b
static const uint64_t VT_DEL            = 4; // 0100b
static const uint64_t VT_BLOCK          = 8; // 1000b -- fourth bit indicates blocksubst
static const uint64_t VT_BLOCKSUBST     = 9; // 1001b
static const uint64_t VT_COMPLEXINS     = 10;// 1010b
static const uint64_t VT_COMPLEXDEL     = 12;// 1100b

/** counts and what they mean */
static const int VS_COUNTS = 256;

/** all alleles / refvars seen */
static const uint64_t CT_ALLELES = 0;

/** locations / input VCF lines */
static const uint64_t CT_LOCATIONS = 1;

/** alleles by type - see above VT_... - 16 different types possible */
static const uint64_t CT_ALLELES_BY_TYPE = 8;

/** alleles by type */
static const uint64_t CT_LOCATIONS_BY_TYPE = 24;

static const char * CT_NAMES [] = {
    "alleles",              // 0
    "locations",
    "homref",
    "haploid",
    "het",                  // 4
    "homalt",
    "hetalt",
    "unknown_gt",
    "unknown_allele",       // 8 + ...
    "snp",                  //      0001
    "ins",                  //      0010
    "snp_ins",              //      0011
    "del",                  //      0100
    "snp_del",              //      0101
    "ins_del",              //      0110
    "snp_ins_del",          //      0111
    "block",                //      1000
    "blocksubst",           //      1001
    "blocksubst_ins",       //      1010
    "blocksubst_snp_ins",   //      1011
    "blocksubst_del",       //      1100
    "blocksubst_snp_del",   //      1101
    "blocksubst_ins_del",   //      1110
    "blocksubst_snp_ins_del"//      1111
    // ... now come combinations of 8+x and 4+y
};


struct VariantStatisticsImpl
{
    VariantStatisticsImpl(const char * ref_fasta, bool _count_homref) : ref(ref_fasta), count_homref(_count_homref),
        alignment(makeAlignment("klibg"))
    {
        memset(counts, 0, sizeof(size_t)*VS_COUNTS);
    }

    ~VariantStatisticsImpl() {}

    VariantStatisticsImpl(VariantStatisticsImpl const & rhs) : ref(rhs.ref),
        alignment(makeAlignment("klibg"))
    {
        memcpy(counts, rhs.counts, sizeof(size_t)*VS_COUNTS);
        count_homref = rhs.count_homref;
    }

    size_t counts[VS_COUNTS];

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
    for (size_t i = 0; i < CT_LOCATIONS_BY_TYPE; ++i)
    {
        if (root.isMember(CT_NAMES[i]))
        {
            _impl->counts[i] = root[CT_NAMES[i]].asUInt64();
        }
        else
        {
            _impl->counts[i] = 0;
        }
    }
    // we have 16 composite allele types, 5 location types.
    for (size_t i = CT_LOCATIONS_BY_TYPE; i < CT_LOCATIONS_BY_TYPE + 16*5; ++i)
    {
        size_t vt = (i - CT_LOCATIONS_BY_TYPE) & 0x0f;
        size_t lt = ((i - CT_LOCATIONS_BY_TYPE) >> 4) & 0x07;
        std::string n(CT_NAMES[CT_ALLELES_BY_TYPE + vt]);
        n += "__";
        n += CT_NAMES[CT_LOCATIONS + 1 + lt];
        if (root.isMember(n.c_str()))
        {
            _impl->counts[i] = root[n.c_str()].asUInt64();
        }
        else
        {
            _impl->counts[i] = 0;
        }
#ifdef DEBUG_VARIANTSTATISTICS
        std::cerr << i << " /  " << n << " : " << _impl->counts[i];
#endif
    }
}

Json::Value VariantStatistics::write()
{
    Json::Value root;
    for (size_t i = 0; i < CT_LOCATIONS_BY_TYPE; ++i)
    {
        if (_impl->counts[i])
        {
            root[CT_NAMES[i]] = Json::Value::UInt64(_impl->counts[i]);
        }
    }
    // we have 16 composite allele types, 5 location types.
    for (size_t i = CT_LOCATIONS_BY_TYPE; i < CT_LOCATIONS_BY_TYPE + 16*5; ++i)
    {
        if (_impl->counts[i])
        {
            size_t vt = (i - CT_LOCATIONS_BY_TYPE) & 0x0f;
            size_t lt = ((i - CT_LOCATIONS_BY_TYPE) >> 4) & 0x07;
            std::string n(CT_NAMES[CT_ALLELES_BY_TYPE + vt]);
            if (std::string(CT_NAMES[CT_LOCATIONS + 1 + lt]) == "homref")
            {
                continue;
            }
            n += "__";
            n += CT_NAMES[CT_LOCATIONS + 1 + lt];
            root[n.c_str()] = Json::Value::UInt64(_impl->counts[i]);
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

void VariantStatistics::add(Variants const & rhs, int sample)
{
    size_t total_als = 0;
    size_t total_snp = 0;
    size_t total_del = 0;
    size_t total_ins = 0;

    auto add_rv = [this, rhs, &total_als, &total_snp, &total_del, &total_ins](RefVar const & rv)
    {
        ++total_als;
        size_t tmp;
        realignRefVar(_impl->ref, rhs.chr.c_str(), rv, _impl->alignment.get(),
                      total_snp, total_ins, total_del, tmp);
    };
    // count hom-alt alleles only once
    if(rhs.calls[sample].isHomalt())
    {
        add_rv(rhs.variation[rhs.calls[sample].gt[0]-1]);
    }
    else
    {
        for(size_t i = 0; i < rhs.calls[sample].ngt; ++i)
        {
            if (rhs.calls[sample].gt[i] > 0)
            {
                add_rv(rhs.variation[rhs.calls[sample].gt[i]-1]);
            }
        }
    }
    for(auto i : rhs.ambiguous_alleles[sample])
    {
        if (i > 0)
        {
            add_rv(rhs.variation[i-1]);
        }
    }

    _impl->counts[CT_ALLELES] += total_als;
    _impl->counts[CT_ALLELES_BY_TYPE + VT_SNP] += total_snp;
    _impl->counts[CT_ALLELES_BY_TYPE + VT_DEL] += total_del;
    _impl->counts[CT_ALLELES_BY_TYPE + VT_INS] += total_ins;

    gttype gtt = getGTType(rhs.calls[sample]);
    uint64_t t = (((uint64_t)gtt) << 4);
    int vts = 0;
    if(total_snp)
    {
        t |= VT_SNP;
        ++vts;
    }
    if(total_del)
    {
        t |= VT_DEL;
        ++vts;
    }
    if(total_ins)
    {
        t |= VT_INS;
        ++vts;
    }
    if(vts > 1)
    {
        t |= VT_BLOCK;
    }

    if (!_impl->count_homref && (CT_LOCATIONS + (t >> 4) + 1 == 2))
    {
        // don't count homref locations
        return;
    }

    ++_impl->counts[CT_LOCATIONS];
    ++_impl->counts[CT_LOCATIONS + (t >> 4) + 1];
    ++_impl->counts[CT_LOCATIONS_BY_TYPE + t];

}

void VariantStatistics::add(const char * chr, RefVar const & rhs)
{
    ++_impl->counts[CT_ALLELES];
    size_t tmp;
    realignRefVar(_impl->ref, chr, rhs, _impl->alignment.get(),
            _impl->counts[CT_ALLELES_BY_TYPE + VT_SNP],
            _impl->counts[CT_ALLELES_BY_TYPE + VT_INS],
            _impl->counts[CT_ALLELES_BY_TYPE + VT_DEL],
            tmp);
}

}
