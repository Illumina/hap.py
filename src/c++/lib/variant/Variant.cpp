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
 *  \brief Variant implementation
 *
 * \file Variant.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <htslib/synced_bcf_reader.h>
#include "VariantImpl.hh"

#include <cmath>
#include <htslib/vcf.h>

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
        if(v.ngt == 0 || v.isNocall())
        {
            o << ".";
        }
        else
        {
            for (size_t i = 0; i < v.ngt; ++i)
            {
                if(i > 0)
                {
                    o << (v.phased ? "|" : "/");
                }
                o << v.gt[i];
            }
        }

        if(!std::isnan(v.qual) && v.qual > 0)
        {
            o << " " << v.qual ;
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

    float Variants::getQual() const
    {
        float qual = 0;
        for(auto const & c : calls)
        {
            qual += c.qual;
        }
        if(calls.empty())
        {
            return 0;
        }
        else
        {
            return qual / calls.size();
        }
    }

    /** interface to set / get INFO values */
    int Variants::getInfoInt(const char * id) const
    {
        if(!infos.isMember(id))
        {
            return bcf_int32_missing;
        }
        try
        {
            return infos[id].asInt();
        }
        catch(std::runtime_error const & )
        {
            return bcf_int32_missing;
        }
    }

    float Variants::getInfoFloat(const char * id) const
    {
        if(!infos.isMember(id))
        {
            return bcfhelpers::missing_float();
        }
        try
        {
            return infos[id].asFloat();
        }
        catch(std::runtime_error const & )
        {
            return bcfhelpers::missing_float();
        }
    }

    std::string Variants::getInfoString(const char * id) const
    {
        if(!infos.isMember(id))
        {
            return "";
        }
        return infos[id].asString();
    }

    bool Variants::getInfoFlag(const char * id) const
    {
        if(!infos.isMember(id))
        {
            return false;
        }
        return infos[id].asBool();
    }

    void Variants::delInfo(const char * id)
    {
        infos.removeMember(id);
    }

    void Variants::setInfo(const char * id, bool flag)
    {
        infos[id] = Json::Value(flag);
    }

    void Variants::setInfo(const char * id, int val)
    {
        infos[id] = Json::Value(val);
    }

    void Variants::setInfo(const char * id, float val)
    {
        infos[id] = Json::Value(val);
    }

    void Variants::setInfo(const char * id, const char * val)
    {
        infos[id] = Json::Value(val);
    }
}
