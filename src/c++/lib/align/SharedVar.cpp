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
 * Find shared variation between REF, ALT1, and ALT2
 *
 * \file SharedVar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "SharedVar.hh"

#include <memory>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include "Error.hh"
#include "MultipleAlignment.hh"

#define DEBUG_SHAREDVAR

namespace haplotypes
{

struct SharedVarImpl
{
    MultipleAlignment maln;
};

SharedVar::SharedVar()
{
    _impl = new SharedVarImpl();
}

SharedVar::~SharedVar()
{
    delete _impl;
}

/**
 * Given the reference and two alts, make a list of variants.
 * 
 * RefVar Flags will indicate whether the variant is present in alt1 (0x01), alt2 (0x02), or both (0x03)
 */
void SharedVar::splitVar(
    int64_t refoffset,
    std::string const & ref,
    std::string const & alt1,
    std::string const & alt2,
    std::list<variant::RefVar> & vars
)
{
    std::string padded_ref, padded_alt1, padded_alt2;
    _impl->maln.align(ref, alt1, alt2, padded_ref, padded_alt1, padded_alt2);

    int64_t refpos = refoffset;

    for (size_t i = 0; i < padded_ref.size(); ++i)
    {
            char r = padded_ref[i], a1 = padded_alt1[i], a2 = padded_alt2[i];
#ifdef _DEBUG
            char tr = toupper(r);
            if(tr != 'A' && tr != 'C' && tr != 'G' && tr != 'T' && tr != 'N' && tr != '-')
            {
                error("Invalid ref allele %c returned near %i", tr, refoffset);
            }
            tr = toupper(a1);
            if(tr != 'A' && tr != 'C' && tr != 'G' && tr != 'T' && tr != 'N' && tr != '-')
            {
                error("Invalid alt1 allele %c returned near %i", tr, refoffset);
            }
            tr = toupper(a2);
            if(tr != 'A' && tr != 'C' && tr != 'G' && tr != 'T' && tr != 'N' && tr != '-')
            {
                error("Invalid alt2 allele %c returned near %i", tr, refoffset);
            }
#endif
            if(a1 == a2)
            {
                // a1 != ref and a2 != ref; shared variant.
                variant::appendToVarList(refpos, r, a1, vars, 0x03);
            }
            else
            {
                variant::appendToVarList(refpos, r, a1, vars, 0x01);
                variant::appendToVarList(refpos, r, a2, vars, 0x02);
            }
            if(r != '-')
            {
                ++refpos;
            }
    }
}

}