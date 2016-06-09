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
 * \brief Homref interval splitter processor
 *
 * \file VariantCallsOnly.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/VariantCallsOnly.hh"
#include "Error.hh"

#include "helpers/IntervalBuffer.hh"
#include "helpers/LocationInfo.hh"

#include <queue>
#include <vector>
#include <list>

/* #define DEBUG_VARIANTCALLSONLY */

namespace variant {

struct VariantCallsOnly::VariantCallsOnlyImpl
{
    VariantCallsOnlyImpl() {}

    VariantCallsOnlyImpl(VariantCallsOnly::VariantCallsOnlyImpl const & rhs)
    {
        buffered_variants = rhs.buffered_variants;
        vs = rhs.vs;
        homref_ivs = rhs.homref_ivs;
        homref_dp = rhs.homref_dp;
    }

    std::list<Variants> buffered_variants;

    // homref information
    intervals::IntervalBuffer homref_ivs;
    std::vector< intervals::LocationInfo<int> > homref_dp;

    Variants vs;
};


VariantCallsOnly::VariantCallsOnly()
{
    _impl = new VariantCallsOnlyImpl();
}

VariantCallsOnly::VariantCallsOnly(VariantCallsOnly const & rhs)
{
    _impl = new VariantCallsOnlyImpl(*rhs._impl);
}

VariantCallsOnly::~VariantCallsOnly()
{
    delete _impl;
}

VariantCallsOnly const & VariantCallsOnly::operator=(VariantCallsOnly const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    delete _impl;
    _impl = new VariantCallsOnlyImpl(*rhs._impl);
    return *this;
}

/** enqueue a set of variants */
void VariantCallsOnly::add(Variants const & v)
{
    if (_impl->buffered_variants.size() > 0 &&
        v.chr == _impl->buffered_variants.back().chr &&
        v.pos < _impl->buffered_variants.back().pos)
    {
        error("Variant added out of order at %s:%i / %i", v.chr.c_str(), v.pos, _impl->buffered_variants.back().pos);
    }
    //
    // keep fails
    if(v.getInfoFlag("IMPORT_FAIL"))
    {
#ifdef DEBUG_VARIANTCALLSONLY
        std::cerr << "fail-pass-on: " << v << "\n";
#endif
        _impl->buffered_variants.push_back(v);
        return;
    }
#ifdef DEBUG_VARIANTCALLSONLY
    std::cerr << "Variants added: " << v << "\n";
#endif
    if (v.anyHomref())
    {
        Variants non_hr = v;
        int n_non_hr = (int) v.calls.size();
        for (size_t q = 0; q < v.calls.size(); ++q)
        {
            if(v.calls[q].isHomref())
            {
                non_hr.calls[q] = Call();
                _impl->homref_ivs.addInterval(v.pos, v.pos + v.len - 1, q);

                // remember dp
                if(q >= _impl->homref_dp.size())
                {
                    _impl->homref_dp.resize(q+1);
                }
                _impl->homref_dp[q].set(v.calls[q].dp, v.pos, v.pos + v.len - 1);
                --n_non_hr;
            }
            else if(v.calls[q].isNocall())
            {
                --n_non_hr;
            }
        }
        if (n_non_hr || non_hr.anyAmbiguous())
        {
#ifdef DEBUG_VARIANTCALLSONLY
            std::cerr << "non-hr-add: " << v << "\n";
#endif
            _impl->buffered_variants.push_back(non_hr);
        }
    }
    else
    {
#ifdef DEBUG_VARIANTCALLSONLY
        std::cerr << "non-hr-pass-on: " << v << "\n";
#endif
        _impl->buffered_variants.push_back(v);
    }
}

/**
 * @brief Return variant block at current position
 **/
Variants & VariantCallsOnly::current()
{
#ifdef DEBUG_VARIANTCALLSONLY
        std::cerr << "VCO output: " << _impl->vs << "\n";
#endif
    return _impl->vs;
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantCallsOnly::advance()
{
    if (_impl->buffered_variants.empty())
    {
        return false;
    }
    _impl->vs = _impl->buffered_variants.front();
    _impl->buffered_variants.pop_front();

    // we return sorted variants, so we can forget homref information before here
    _impl->homref_ivs.advance(_impl->vs.pos-1);
    for (auto & dp : _impl->homref_dp)
    {
        dp.reset_to(_impl->vs.pos-1);
    }

    // fix homref calls if output variant is covered
    for (size_t q = 0; q < _impl->vs.calls.size(); ++q) {
        if(_impl->vs.calls[q].isNocall() &&
           _impl->homref_ivs.isCovered(_impl->vs.pos, _impl->vs.pos + _impl->vs.len - 1, q))
        {
            _impl->vs.calls[q].ngt = 2;
            _impl->vs.calls[q].gt[0] = 0;
            _impl->vs.calls[q].gt[1] = 0;
            if (_impl->homref_dp.size() > q)
            {
                _impl->homref_dp[q].queryMean(_impl->vs.pos,
                                              _impl->vs.pos + _impl->vs.len - 1,
                                              _impl->vs.calls[q].dp);
            }
        }
    }

    return true;
}

/** empty internal buffer */
void VariantCallsOnly::flush()
{
    _impl->buffered_variants.clear();
    _impl->vs = Variants();
    _impl->homref_ivs.advance(-1);
}

}
