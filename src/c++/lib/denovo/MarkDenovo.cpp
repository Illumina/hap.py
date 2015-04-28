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
 * \brief 
 *
 * \file MarkDenovo.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#include "denovo/MarkDenovo.hh"

#include "IntervalTree.h"

#include <set>

namespace variant
{

struct MarkDenovo::MarkDenovoImpl
{
    MarkDenovoImpl() : firstone(false), added_new(false) {}
    MarkDenovoImpl(MarkDenovoImpl const & rhs) : buffer(rhs.buffer), tmp(rhs.tmp), firstone(rhs.firstone) {}
    std::list<Variants> buffer; 
    Variants tmp;
    bool firstone;
    bool added_new;
};

MarkDenovo::MarkDenovo()
{
    _impl = new MarkDenovoImpl();
}

MarkDenovo::~MarkDenovo()
{
    delete _impl;
}

MarkDenovo::MarkDenovo(MarkDenovo const & rhs)
{
    _impl = new MarkDenovoImpl(*rhs._impl);
}

MarkDenovo const & MarkDenovo::operator=(MarkDenovo const & rhs)
{
    if(&rhs != this)
    {
        delete _impl;
        _impl = new MarkDenovoImpl(*rhs._impl);
    }
    return *this;
}


/** Variant input **/
/** enqueue a set of variants */
void MarkDenovo::add(Variants const & vs) 
{ 
    _impl->added_new = true;
    if (_impl->buffer.empty())
    {
        _impl->firstone = true;
    }
    _impl->buffer.push_back(vs); 
}

/** Variant output **/
/**
 * @brief Return variant block at current position
 **/
Variants & MarkDenovo::current() 
{
    if( _impl->buffer.empty() ) 
    { 
        return _impl->tmp; 
    } 
    else 
    {
        return _impl->buffer.front(); 
    }  
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool MarkDenovo::advance() 
{
    if(_impl->added_new  && !_impl->buffer.empty())
    {
        _impl->added_new = false;

        std::string chr;
        int64_t start = -1, end = -1;

        int n_nonsnp = 0;

        for (Variants & v : _impl->buffer)
        {
            if(chr.size() == 0)
            {
                chr = v.chr;
            }
            if(chr != v.chr)
            {
                error("Error in MarkDenovo block processing: Use VariantProcessor in blocked mode too prevent issues at chr boundary.");
            }

            if(v.calls.size() < 3)
            {
                error("Error in MarkDenovo block processing: I need at least three samples.");
            }

            if(start < 0)
            {
                start = v.pos;
            }
            else
            {
                start = std::min(start, v.pos);
            }

            if(end < 0)
            {
                end = v.pos + v.len - 1;
            }
            else
            {
                end = std::max(end, v.pos + v.len - 1);
            }

            bool ns_found = false;
            for (Call & c : v.calls)
            {
                for (size_t i = 0; i < c.ngt; ++i)
                {
                    if(c.gt[i] > 0)
                    {
                        int64_t reflen = v.variation[c.gt[i]-1].end - v.variation[c.gt[i]-1].start + 1;
                        int64_t altlen = v.variation[c.gt[i]-1].alt.size();
                        if(reflen != 1 || altlen != 1)
                        {
                            ns_found = true;
                            break;
                        }
                    }
                }
                if (ns_found)
                {
                    ++n_nonsnp;
                    break;
                }
            }
        }

        // if(n_nonsnp)
        // {
        //     for (Variants & v : _impl->buffer)
        //     {
        //         setVariantInfo(v, "denovo", "");
        //     }
        // }
        // else
        // {
            // Denovo SNPs are easy, we do these location by location
            for (Variants & v : _impl->buffer)
            {
                Call & child = v.calls[0];
                Call & mum  = v.calls[1];
                Call & dad  = v.calls[2];

                if (mum.ngt == 0 || dad.ngt == 0 || child.ngt == 0)
                {
                    setVariantInfo(v, "denovo", "");
                    continue;                    
                }

                bool any_nocall = false;
                std::set<int> parental_gts;
                for (size_t j = 0; j < mum.ngt; ++j)
                {
                    if(mum.gt[j] >= 0)
                    {
                        parental_gts.insert(mum.gt[j]);
                    }
                    else
                    {
                        any_nocall = true;
                    }
                }

                if(any_nocall)
                {
                    setVariantInfo(v, "denovo", "");
                    continue;
                }

                for (size_t j = 0; j < dad.ngt; ++j)
                {
                    if(dad.gt[j] >= 0)
                    {
                        parental_gts.insert(dad.gt[j]);
                    }
                    else
                    {
                        any_nocall = true;
                    }
                }

                if(any_nocall)
                {
                    setVariantInfo(v, "denovo", "");
                    continue;
                }

                bool parents_homref = parental_gts.size() == 1 && parental_gts.count(0);

                for(size_t i = 3; i < v.calls.size(); ++i)
                {
                    Call & sib = v.calls[i];
                    for (size_t j = 0; j < sib.ngt; ++j)
                    {
                        if(sib.gt[j] >= 0)
                        {
                            parental_gts.insert(sib.gt[j]);
                        }
                    }
                }

                bool dn_found = false;
                bool child_only_homref = true;

                for (size_t j = 0; j < child.ngt; ++j)
                {
                    if(child.gt[j] == 0)
                    {
                        if(!parental_gts.count(0))
                        {
                            dn_found = true;
                        }
                    }
                    else if(child.gt[j] > 0)
                    {
                        child_only_homref = false;
                        if(!parental_gts.count(child.gt[j]))
                        {
                            dn_found = true;
                        }
                    }
                    else
                    {
                        any_nocall = true;
                    }
                }                

                if(any_nocall)
                {
                    setVariantInfo(v, "denovo", "");
                    continue;
                }

                if(dn_found)
                {
                    if(parents_homref)
                    {
                        setVariantInfo(v, "denovo", "child_only");
                    }
                    else
                    {
                        setVariantInfo(v, "denovo", child_only_homref ? "parents_only" : "child_partial");
                    }
                }
            }
        // }
    }

    if (_impl->firstone && !_impl->buffer.empty())
    {
        _impl->firstone = false;
        return true;
    }
    else
    {
        if(!_impl->buffer.empty()) 
        {
            _impl->buffer.pop_front();
        }
        return !_impl->buffer.empty(); 
    }
}

/** empty internal buffer */
void MarkDenovo::flush() 
{ 
    _impl->buffer.clear(); 
}

}
