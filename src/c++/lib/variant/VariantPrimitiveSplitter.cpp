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
 * \brief Class to split complex alleles into primitive blocks
 *
 * \file VariantPrimitiveSplitter.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantPrimitiveSplitter.hh"

#include <iostream>
#include <queue>
#include <vector>
#include <memory>

#include "Alignment.hh"

/* #define DEBUG_VARIANTPRIMITIVESPLITTER */

namespace variant {


namespace ns_variantprimitivesplitter
{
    struct VariantCompare
    {
        bool operator() (Variants const & v1, Variants const & v2)
        {
            return v1.pos > v2.pos;
        }
    };

    typedef std::priority_queue<
        Variants,
        std::vector<Variants>,
        ns_variantprimitivesplitter::VariantCompare
    > VariantQueue;
}

struct VariantPrimitiveSplitter::VariantPrimitiveSplitterImpl
{
    VariantPrimitiveSplitterImpl() : aln(makeAlignment("klibg")) {}

    VariantPrimitiveSplitterImpl(VariantPrimitiveSplitter::VariantPrimitiveSplitterImpl const & rhs)
    	: buffered_variants(rhs.buffered_variants), output_variants (rhs.output_variants), vs(rhs.vs),
          reference(rhs.reference), aln(makeAlignment("klibg"))
    {
        if (reference != "")
        {
            ref_fasta = std::move( std::unique_ptr<FastaFile>(new FastaFile( reference.c_str() ) ) );
        }
    }

    std::vector<Variants> buffered_variants;
    ns_variantprimitivesplitter::VariantQueue output_variants;

    Variants vs;

    std::string reference;
    std::unique_ptr<FastaFile> ref_fasta;
    std::unique_ptr<Alignment> aln;

    /** Push a variant to output_variants
     *
     *  This is somewhat comples because we need to take apart a complex variant
     *  into alt and homref calls here.
     *  Some of the alleles might be deletions, we return the homref bits
     *  independently.
     */
    void pushOutputVariant(Variants vs)
    {
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
        std::cerr << "pov in: " << vs << "\n";
#endif
        for(size_t ci = 0; ci < vs.calls.size(); ++ci)
        {
            Call & c = vs.calls[ci];
            for(size_t gti = 0; gti < c.ngt; ++gti)
            {
                if(c.gt[gti] > 0)
                {
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                    assert(vs.pos == std::min(vs.pos, vs.variation[c.gt[gti]-1].start));
#endif
                    vs.len = std::max(vs.len, vs.variation[c.gt[gti]-1].end - vs.variation[c.gt[gti]-1].start + 1);
                }
            }
        }
        // if length is 1, we can report the homref calls together with the SNPs
        size_t homref_calls_total = 0, calls_total = 0;
        if(vs.anyHomref())
        {
            if(vs.len == 1)
            {
                for(size_t ci = 0; ci < vs.calls.size(); ++ci)
                {
                    if(!vs.calls[ci].isNocall())
                    {
                        ++calls_total;
                    }
                    if(vs.calls[ci].isHomref() && vs.calls[ci].isHemi())
                    {
                        vs.calls[ci] = Call();
                        homref_calls_total++;
                    }
                }
                if(homref_calls_total != calls_total)
                {
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                    std::cerr << "out-push: " << vs << "\n";
#endif
                    output_variants.push(vs);
                }
            }
            else
            {
                Variants vs_homref = vs;
                int calls_added = 0;
                for(size_t ci = 0; ci < vs.calls.size(); ++ci)
                {
                    if(!vs.calls[ci].isNocall())
                    {
                        ++calls_total;
                    }
                    if(vs.calls[ci].isHomref())
                    {
                        if(vs.calls[ci].isHemi())
                        {
                            vs_homref.calls[ci] = Call();
                        }
                        else
                        {
                            ++calls_added;
                        }
                        vs.calls[ci] = Call();
                        ++homref_calls_total;
                    }
                    else
                    {
                        vs_homref.calls[ci] = Call();
                    }
                }
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                std::cerr << "homref: " << vs_homref << "\n";
#endif
                if(calls_added > 0)
                {
                    vs_homref.len = 1;
                    vs_homref.variation.clear();
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                    std::cerr << "out-push: " << vs_homref << "\n";
#endif
                    output_variants.push(vs_homref);
                }
                // only add if we haven't actually added everything above already
                if(homref_calls_total != calls_total)
                {
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                    std::cerr << "out-push: " << vs << "\n";
#endif
                    output_variants.push(vs);
                }
            }
        }
        else
        {
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
            std::cerr << "out-push: " << vs << "\n";
#endif
            output_variants.push(vs);
        }
    }
};



VariantPrimitiveSplitter::VariantPrimitiveSplitter() : _impl(new VariantPrimitiveSplitterImpl())
{
}

VariantPrimitiveSplitter::VariantPrimitiveSplitter(VariantPrimitiveSplitter const & rhs) : _impl(
    new VariantPrimitiveSplitterImpl(*rhs._impl))
{
}

VariantPrimitiveSplitter::~VariantPrimitiveSplitter()
{
    delete _impl;
}

VariantPrimitiveSplitter & VariantPrimitiveSplitter::operator=(VariantPrimitiveSplitter const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    delete _impl;
    _impl = new VariantPrimitiveSplitterImpl(*rhs._impl);
    return *this;
}

/**
 * @brief Reference for shifting
 *
 */
void VariantPrimitiveSplitter::setReference(std::string const & fasta)
{
    _impl->reference = fasta;
    _impl->ref_fasta = std::move( std::unique_ptr<FastaFile>(new FastaFile(_impl->reference.c_str() ) ) );
}


/** enqueue a set of variants */
void VariantPrimitiveSplitter::add(Variants const & vs)
{
    if (_impl->buffered_variants.size() > 0 &&
        vs.chr == _impl->buffered_variants.back().chr &&
        vs.pos < _impl->buffered_variants.back().pos)
    {
        error("Variant added out of order at %s:%i / %i", vs.chr.c_str(), vs.pos, _impl->vs.pos);
    }
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
    std::cerr << "VHPS input: " << vs << "\n";
#endif
    _impl->buffered_variants.push_back(vs);
}

/**
 * @brief Return variant block at current position
 **/
Variants & VariantPrimitiveSplitter::current()
{
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
    for(auto const & c : _impl->vs.calls)
    {
        for(size_t j = 0; j < c.ngt; ++j)
        {
            assert(c.gt[j] <= (signed)_impl->vs.variation.size());
        }
    }
    std::cerr << "VHPS-output: " << _impl->vs << "\n";
#endif
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantPrimitiveSplitter::advance()
{
    // refill output buffer?
    if (_impl->output_variants.empty())
    {
        for (Variants & v : _impl->buffered_variants)
        {
            // homref?
            if(v.variation.size() == 0) {
                _impl->output_variants.push(v);
                continue;
            }
            bool any_realignable = false;
            for (size_t i = 0; i < v.variation.size(); ++i)
            {
                int64_t reflen = v.variation[i].end - v.variation[i].start + 1;
                int64_t altlen = (int64_t)v.variation[i].alt.size();
                if (reflen > 1 && altlen > 0)
                {
                    any_realignable = true;
                    break;
                }
            }

            if(!any_realignable)
            {
                _impl->output_variants.push(v);
                continue;
            }

            // remove all homref calls and add them
            Variants v_homref = v;
            // these get passed on by v_homref
            for(auto & l : v.ambiguous_alleles)
            {
                l.clear();
            }

            bool any_homref = false;
            for (size_t i = 0; i < v.calls.size(); ++i)
            {
                if (v.calls[i].isHomref())
                {
                    any_homref = true;
                    v.calls[i] = Call();
                }
                else
                {
                    v_homref.calls[i] = Call();
                }
            }

            // push homref / ambiguous variant parts
            if (any_homref || v_homref.anyAmbiguous())
            {
                _impl->output_variants.push(v_homref);
            }

            ns_variantprimitivesplitter::VariantQueue output_queue;

            std::vector< std::list<RefVar> > rvlists(v.variation.size());
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
            std::cerr << "at " << v.pos  << "\n";
#endif
            for (size_t i = 0; i < v.variation.size(); ++i)
            {
                realignRefVar(*(_impl->ref_fasta), v.chr.c_str(), v.variation[i],
                              _impl->aln.get(), rvlists[i]);
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                for(auto const & rv : rvlists[i])
                {
                    std::cerr << "realigned RV " << i << " : " << rv << "\n";
                }
#endif
            }

            // this really splits calls along samples as well as along
            // primitives. We run a simplified merge-by-location after this
            // to fix things up such that the same calls line up again
            Variants vnew;
            vnew.chr = v.chr;
            vnew.info = v.info;
            vnew.calls.resize(v.calls.size());
            for(size_t ci = 0; ci < v.calls.size(); ++ci)
            {
                if(ci > 0)
                {
                    v.calls[ci-1] = Call();
                }
                Call & c = v.calls[ci];

                if (c.isNocall())
                {
                    continue;
                }

                for(size_t xci = 0; xci < vnew.calls.size(); ++xci)
                {
                    vnew.calls[xci] = Call();
                }

                for(size_t gti = 0; gti < c.ngt; ++gti)
                {
                    int rvallele = -1;
                    int dp = 0;
                    if(c.gt[gti] < 0)
                    {
                        continue;
                    }
                    else if(c.gt[gti] == 0)
                    {
                        dp = c.ad_ref;
                    }
                    else if(c.gt[gti] > 0)
                    {
                        rvallele = c.gt[gti] - 1;
                        dp = c.ad[gti];
                    }

                    Call refCall;
                    refCall.ngt = 1;
                    refCall.gt[0] = 0;
                    refCall.ad[0] = dp;
                    refCall.ad_ref = c.ad_ref;
                    refCall.ad_other = c.ad_other;
                    for(size_t ff = 0; ff < c.nfilter; ++ff)
                    {
                        refCall.filter[ff] = c.filter[ff];
                    }
                    refCall.nfilter = c.nfilter;

                    Call altCall;
                    altCall.ngt = 1;
                    altCall.gt[0] = 1;
                    altCall.ad[0] = dp;
                    altCall.ad_ref = c.ad_ref;
                    altCall.ad_other = c.ad_other;
                    for(size_t ff = 0; ff < c.nfilter; ++ff)
                    {
                        altCall.filter[ff] = c.filter[ff];
                    }
                    altCall.nfilter = c.nfilter;

                    std::list<RefVar>::const_iterator rvx;
                    if((size_t)rvallele < rvlists.size())
                    {
                        rvx = rvlists[rvallele].begin();
                    }
                    int64_t start = v.pos, end = v.pos + v.len - 1;
                    while(start <= end)
                    {
                        if((size_t)rvallele < rvlists.size() &&
                            rvx != rvlists[rvallele].end() && start >= std::min(rvx->end, rvx->start))
                        {
                            vnew.variation.resize(1);
                            vnew.variation[0] = *rvx;
                            vnew.pos = rvx->start;
                            vnew.len = std::max(int64_t(1), rvx->end - rvx->start + 1);
                            vnew.calls[ci] = altCall;
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                            std::cerr << "pushing : " << vnew << "\n";
#endif
                            output_queue.push(vnew);
                            start = rvx->end + 1;
                            ++rvx;
                        }
                        else
                        {
                            vnew.pos = start;
                            vnew.len = 1;
                            vnew.len = 1;
                            vnew.variation.clear();
                            vnew.calls[ci] = refCall;
#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                            std::cerr << "pushing : " << vnew << "\n";
#endif
                            output_queue.push(vnew);

                            start += 1;
                        }
                    }
                }
            }

            Variants vs;
            vs.pos = -1;
            vs.calls.resize(v.calls.size());
            bool has_vs = false;
            // go through output queue and merge calls by location
            //

            while(!output_queue.empty())
            {
                Variants vs2 = output_queue.top();
                output_queue.pop();
                if(vs2.pos != vs.pos || !has_vs)
                {
                    if(has_vs)
                    {
                        _impl->pushOutputVariant(vs);
                    }
                    vs = vs2;
                    has_vs = true;
                    continue;
                }

#ifdef DEBUG_VARIANTPRIMITIVESPLITTER
                assert(vs2.calls.size() == vs.calls.size());
                assert(vs2.calls.size() == v.calls.size());
#endif
                size_t gt_ofs = vs.variation.size();
                for(auto const & var : vs2.variation)
                {
                    vs.variation.push_back(var);
                }

                // go through calls. If the positions match, we should be able to merge
                for(size_t ci = 0; ci < vs.calls.size(); ++ci)
                {
                    Call & c1 = vs.calls[ci];
                    Call const & c2 = vs2.calls[ci];

                    static const auto combine_filters = [](Call & c1, Call const& c2)
                    {
                        for (size_t q = 0; q < c2.nfilter; ++q)
                        {
                            bool found = false;
                            for (size_t q2 = 0; q2 < c1.nfilter; ++q2)
                            {
                                if(c1.filter[q2] == c1.filter[q])
                                {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found && c1.nfilter < MAX_FILTER)
                            {
                                c1.filter[c1.nfilter] = c2.filter[q];
                                ++c1.nfilter;
                            }
                        }
                    };

                    // no-calls can be overwritten
                    if(c1.isNocall())
                    {
                        c1 = c2;
                        // remap GT
                        for(size_t j = 0; j < c1.ngt; ++j)
                        {
                            if(c1.gt[j] > 0)
                            {
                                c1.gt[j] += gt_ofs;
                            }
                        }
                        continue;
                    }

                    if(c2.isNocall())
                    {
                        // ignore nocalls
                        continue;
                    }

                    // merge two hemi-alleles into one het(-alt)
                    if (c1.isHemi() && c2.isHemi())
                    {
                        c1.ngt = 2;
                        c1.gt[1] = c2.gt[0];
                        if(c1.gt[1] > 0)
                        {
                            c1.gt[1] += gt_ofs;
                        }
                        c1.ad[1] = c2.ad[0];
                        combine_filters(c1, c2);
                        continue;
                    }

                    // now that we know that c2 isn't a no-call, we can overwrite c1 with it if c1 is homref
                    if(c1.isHomref())
                    {
                        c1 = c2;
                        // remap GT
                        for(size_t j = 0; j < c1.ngt; ++j)
                        {
                            if(c1.gt[j] > 0)
                            {
                                c1.gt[j] += gt_ofs;
                            }
                        }
                        combine_filters(c1, c2);
                        continue;
                    }

                    // c1 isn't no-call or homref, so it will still have the reference information in it. No need to add hr again
                    // except in the case where both are hemi, which is handled above
                    if(c2.isHomref())
                    {
                        combine_filters(c1, c2);
                        continue;
                    }

                    int vgts[MAX_GT*2];
                    int vdps[MAX_GT*2];

                    size_t nvgt = 0, g2 = 0;
                    bool has_ref = false;
                    int dp_ref = c1.ad_ref;
                    while(g2 < c1.ngt)
                    {
                        if(c1.gt[g2] > 0)
                        {
                            vgts[nvgt] = c1.gt[g2];
                            vdps[nvgt] = c1.ad[g2];
                            nvgt++;
                        }
                        else
                        {
                            has_ref = true;
                            dp_ref = c1.ad_ref;
                        }
                        ++g2;
                    }
                    g2 = 0;
                    while(g2 < c2.ngt)
                    {
                        if(c2.gt[g2] > 0)
                        {
                            vgts[nvgt] = c2.gt[g2] + gt_ofs;
                            vdps[nvgt] = c2.ad[g2];
                            nvgt++;
                        }
                        else
                        {
                            has_ref = true;
                            dp_ref = c2.ad_ref;
                        }
                        ++g2;
                    }
                    c1.ad_ref = dp_ref;
                    if(nvgt == 1)
                    {
                        if(has_ref)
                        {
                            c1.ad[0] = dp_ref;
                            c1.gt[0] = 0;
                            c1.ad[1] = vdps[0];
                            c1.gt[1] = vgts[0];
                        }
                        else
                        {
                            c1.ad[0] = vdps[0];
                            c1.gt[0] = 1;
                        }
                    }
                    else if(nvgt == 2)
                    {
                        c1.ad[0] = vdps[0];
                        c1.gt[0] = vgts[0];
                        c1.ad[1] = vdps[1];
                        c1.gt[1] = vgts[1];
                    }
                    else
                    {
                        error("Allele split resolved to ambiguous calls at %s:%i", vs.chr.c_str(), vs.pos);
                    }
                    combine_filters(c1, c2);
                }
            }
            if(has_vs)
            {
                _impl->pushOutputVariant(vs);
            }
        }
        _impl->buffered_variants.clear();
    }

    if (_impl->output_variants.empty())
    {
        return false;
    }
    _impl->vs = _impl->output_variants.top();
    _impl->output_variants.pop();
    return true;
}

/** empty internal buffer */
void VariantPrimitiveSplitter::flush()
{
    _impl->buffered_variants.clear();
    _impl->output_variants = ns_variantprimitivesplitter::VariantQueue();
    _impl->vs = Variants();
}

}
