// -*- mode: c++; indent-tabs-mode: nil; -*-
//
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
 *  \brief Helper functions for getting variants from VCF/BCF
 *
 *
 * \file Variant.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <string>
#include <sstream>
#include <map>
#include <set>

#include "helpers/StringUtil.hh"
#include "helpers/BCFHelpers.hh"

#include "RefVar.hh"

#include "json/json.h"

namespace variant {

#ifndef MAX_FILTER
#define MAX_FILTER 20
#endif

// GT types
enum gttype
{
    gt_homref = 0,
    gt_haploid = 1,
    gt_het = 2,
    gt_homalt = 3,
    gt_hetalt = 4,
    gt_unknown = 5
};

/**
 * @brief Variant call for a given location
 */
struct Call {
    Call() : ad_ref(-1), ad_other(-1), ngt(0), phased(false), nfilter(0), dp(-1), qual(0)
    {
        for(int i = 0; i < MAX_GT; ++i)
        {
            gt[i] = -1;
            ad[i] = -1;
        }
    }

    inline bool isNocall() const
    {
        for (size_t i = 0; i < ngt; ++i)
        {
            if (gt[i] >= 0)
            {
                return false;
            }
        }
        return true;
    }

    inline bool isHalfcall() const
    {
        return ngt == 2 && ( (gt[0] >= 0 && gt[1] < 0) || (gt[1] >= 0 && gt[0] < 0));
    }

    inline bool isHomref() const
    {
        for (size_t i = 0; i < ngt; ++i)
        {
            if (gt[i] != 0)
            {
                return false;
            }
        }
        return ngt > 0;
    }

    inline bool isHet() const
    {
        return ngt == 2 && ((gt[0] == 0 && gt[1] > 0) || (gt[0] > 0 && gt[1] == 0));
    }

    inline bool isHetAlt() const
    {
        return ngt == 2 && (gt[0] > 0) && (gt[1] > 0) && (gt[0] != gt[1]);
    }

    inline bool isHomalt() const
    {
        return ngt == 2 && gt[0] == gt[1] && gt[1] > 0;
    }

    inline bool isHemi() const
    {
        return ngt == 1;
    }

    int gt[MAX_GT];

    // we keep ref and sum of other ADs + AD for the alleles we see
    int ad_ref, ad_other;
    int ad[MAX_GT];

    size_t ngt;
    bool phased;

    std::string filter[MAX_FILTER];
    size_t nfilter;

    int dp;

    float qual;

    Json::Value formats;
};

/**
 * @brief Classify a variant's GT type
 *
 */
extern gttype getGTType(Call const& var);

/**
 * @brief Stores multiple VCF/BCF Variant records
 * @details Performs basic validation
 *
 */
struct Variants
{
    Variants();

    // variant ordering by creation time
    uint64_t id;
    static uint64_t MAX_VID;

    std::string chr;

    std::vector<RefVar> variation;

    // these are the resulting variant calls
    // for the given location
    std::vector<Call> calls;

    // these capture the extent of the alleles in variation
    int64_t pos;
    int64_t len;

    // We collect all the alleles called in each sample.
    // This captures cases where cannot resolve a diploid
    // genotype
    std::vector< std::list<int> > ambiguous_alleles;

    // Store INFO entries
    Json::Value infos;

    /* return if any calls are homref */
    inline bool anyHomref() const {
        for(Call const & c : calls) {
            if (c.isHomref())
            {
                return true;
            }
        }
        return false;
    }

    /* return if all calls are homref */
    inline bool allHomref() const {
        for(Call const & c : calls) {
            if (!c.isHomref())
            {
                return false;
            }
        }
        return calls.size() > 0;
    }

    /* return if any calls are homref */
    inline bool anyAmbiguous() const {
        for(std::list<int> const & c : ambiguous_alleles) {
            if (!c.empty())
            {
                return true;
            }
        }
        return false;
    }

    float getQual() const;

    /** interface to set / get INFO values */
    int getInfoInt(const char * id) const;
    float getInfoFloat(const char * id) const;
    std::string getInfoString(const char * id) const;
    bool getInfoFlag(const char * id) const;

    void delInfo(const char * id);
    void setInfo(const char * id, bool);
    void setInfo(const char * id, int);
    void setInfo(const char * id, float);
    void setInfo(const char * id, const char *);
};

/** variant comparison operator that makes sure indels go after SNPs
 *  and that insertions go last (because they get ref-padded and actually
 *  start after everything else) */
struct VariantCompare
{
    bool operator() (Variants const * v1, Variants const * v2)
    {
        return (*this)(*v1, *v2);
    }

    bool operator() (Variants const & v1, Variants const & v2)
    {
        if(v1.pos == v2.pos)
        {
            // make sure indels come after SNPs
            int v1_types = 0;
            int v2_types = 0;
            size_t v1_max_reflen = 0;
            size_t v2_max_reflen = 0;
            size_t v1_max_altlen = 0;
            size_t v2_max_altlen = 0;
            for(RefVar const & rv : v1.variation)
            {
                int64_t reflen = rv.end - rv.start + 1;
                v1_max_reflen = std::max(reflen, (int64_t)v1_max_reflen);
                v1_max_altlen = std::max(rv.alt.size(), v1_max_altlen);
                if(rv.end - rv.start == 0 && rv.alt.size() == 1)
                {
                    v1_types = v1_types | 4;
                }
                else if(reflen == 0 && rv.alt.size() > 0)
                {
                    // insertions go last
                    v1_types = v1_types | 1;
                }
                else // non-ref-padded stuff goes after snps, but before ref-padded variants
                {
                    v1_types = v1_types | 2;
                }
            }
            for(RefVar const & rv : v2.variation)
            {
                int64_t reflen = rv.end - rv.start + 1;
                v2_max_reflen = std::max(rv.end - rv.start + 1, (int64_t)v2_max_reflen);
                v2_max_altlen = std::max(rv.alt.size(), v2_max_altlen);
                if(rv.end - rv.start == 0 && rv.alt.size() == 1)
                {
                    v2_types = v2_types | 4;
                }
                else if(reflen == 0 && rv.alt.size() > 0)
                {
                    // insertions go last
                    v2_types = v2_types | 1;
                }
                else // non-ref-padded stuff goes after snps, but before ref-padded variants
                {
                    v2_types = v2_types | 2;
                }
            }

            if(v1_types == v2_types)
            {
                if(v1_max_reflen == v2_max_reflen)
                {
                    if (v1_max_altlen == v2_max_altlen)
                    {
                        std::vector<std::string> v1_alts;
                        std::vector<std::string> v2_alts;
                        // lexicographic sort on alts
                        for(RefVar const & rv : v1.variation)
                        {
                            v1_alts.push_back(rv.alt);
                        }
                        for(RefVar const & rv : v2.variation)
                        {
                            v2_alts.push_back(rv.alt);
                        }
                        std::sort(v1_alts.begin(), v1_alts.end());
                        std::sort(v2_alts.begin(), v2_alts.end());
                        std::string v1a;
                        for (auto s : v1_alts) {
                            v1a += s;
                            v1a += ";";
                        }
                        std::string v2a;
                        for (auto s : v2_alts) {
                            v2a += s;
                            v2a += ";";
                        }

                        if(v1a == v2a)
                        {
                            return v1.id >= v2.id;
                        }
                        else
                        {
                            return v1a >= v2a;
                        }
                    }
                    else
                    {
                        return v1_max_altlen >= v2_max_altlen;
                    }
                }
                else
                {
                    return v1_max_reflen > v2_max_reflen;
                }
            }
            else
            {
                return v1_types < v2_types;
            }
        }
        else
        {
            return v1.pos > v2.pos;
        }
    }
};

extern std::ostream & operator<<(std::ostream &o, gttype const & v);
extern std::ostream & operator<<(std::ostream &o, Call const & v);
extern std::ostream & operator<<(std::ostream &o, Variants const & v);

enum class VariantBufferMode {
    // buffer single variant records
    buffer_count,
    // buffer up to block overlap boundary -- param is the
    // number of BP we need to be clear of the last RefVar end
    // to start a new block
    buffer_block,
    // read entire file into memory -- param isn't used
    buffer_all,
    // buffer up to position
    buffer_endpos,
};

/** implementation is in VariantProcessor.cpp */
class VariantReader;
class AbstractVariantProcessingStep
{
public:
    /** Interface **/

    /** Variant input **/
    /** enqueue a set of variants */
    virtual void add(Variants const & vs) = 0;

    /** Variant output **/

    /**
     * @brief Return variant block at current position
     **/
    virtual Variants & current() = 0;

    /**
     * @brief Advance one line
     * @return true if a variant was retrieved, false otherwise
     */
    virtual bool advance() = 0;

    /** Flush internal buffers */
    virtual void flush() = 0;

    /**
     * @brief Add variants from a VariantReader (see above for VariantBufferMode and param)
     */
    bool add(VariantReader & source,
            VariantBufferMode bt,
            int64_t param = 0
           );

    /** add refvar record to a given variant */
    void add_variant(int sample, const char * chr, RefVar const & rv, bool het);

    /** enqueue homref block */
    void add_homref(int sample, const char * chr, int64_t start, int64_t end, bool het);

};


} // namespace variant

#include "variant/VariantReader.hh"
#include "variant/VariantWriter.hh"
#include "variant/VariantLocationMap.hh"
#include "variant/VariantStatistics.hh"

#include "variant/VariantProcessor.hh"
