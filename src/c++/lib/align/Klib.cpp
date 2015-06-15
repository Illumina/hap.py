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
 * \brief klib Smith Waterman wrapped. See external/klib.
 *
 * \file Klib.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Klib.hh"
#include "KlibImpl.hh"

KlibAlignment::KlibAlignment()
{
    _impl = new KlibAlignmentImpl();
    _impl->valid_result = false;
    setParameters(AlignmentParameters());
}

KlibAlignment::~KlibAlignment()
{
    delete _impl;
}


void KlibAlignment::setParameters(AlignmentParameters const & ap)
{
    for (int i = 0; i < 25; ++i)
    {
        _impl->mat[i] = ap.subs_mat[i];
    }
    _impl->gapo = ap.gapo;
    _impl->gape = ap.gape;
}

void KlibAlignment::getParameters(AlignmentParameters & ap)
{
    for (int i = 0; i < 25; ++i)
    {
        ap.subs_mat[i] = _impl->mat[i];
    }
    ap.gapo = _impl->gapo;
    ap.gape = _impl->gape;
}

/**
 * @brief set target sequence
 */
void KlibAlignment::setRef(const char * seq)
{
    _impl->valid_result = false;
    _impl->reflen = strlen(seq);
    _impl->ref = std::shared_ptr<uint8_t>(new uint8_t[_impl->reflen], []( uint8_t *p ) { delete[] p; });
    translate(seq, _impl->ref.get(), _impl->reflen);
}

/*
 * @brief set query sequence
 */
void KlibAlignment::setQuery(const char * seq)
{
    if(_impl->qprofile)
    {
        free(_impl->qprofile);
        _impl->qprofile = NULL;
    }
    _impl->valid_result = false;
    _impl->altlen = strlen(seq);
    _impl->alt = std::shared_ptr<uint8_t>(new uint8_t[_impl->altlen], []( uint8_t *p ) { delete[] p; });
    translate(seq, _impl->alt.get(), _impl->altlen);
}

/**
 * @brief Get the alignment score
 */
int KlibAlignment::getScore()
{
    if(!_impl->valid_result)
    {
        this->update();
    }
    return _impl->result.score;
}

/**
 * @brief Get a cigar string
 */
void KlibAlignment::getCigar(int &r0, int & r1, int &a0, int & a1, std::string &cig)
{
    if(!_impl->valid_result)
    {
        this->update();
    }
    r0 = _impl->result.qb;
    r1 = _impl->result.qe;
    a0 = _impl->result.tb;
    a1 = _impl->result.te;
    cig = makeCigar(_impl->result.tb, _impl->result.te, _impl->altlen, _impl->cigar_len, _impl->cigar);
}

/**
 * @brief Get int* cigar string + start + end
 */
void  KlibAlignment::getCigar(int &r0, int & r1, int &a0, int & a1, int & n_cigar, uint32_t *& cigar)
{
    if(!_impl->valid_result)
    {
        this->update();
    }
    r0 = _impl->result.qb;
    r1 = _impl->result.qe;
    a0 = _impl->result.tb;
    a1 = _impl->result.te;
    n_cigar = _impl->cigar_len;
    cigar = _impl->cigar;
}


/**
 * @brief Debug dump, optional
 */
void KlibAlignment::update()
{
    _impl->result = ksw_align(
        _impl->reflen, _impl->ref.get(),
        _impl->altlen, _impl->alt.get(),
        5, _impl->mat, _impl->gapo, _impl->gape,
        KSW_XSTART,     // add flags here
        &(_impl->qprofile));

    if(_impl->cigar)
    {
        free(_impl->cigar);
        _impl->cigar = NULL;
        _impl->cigar_len = 0;
    }

    ksw_global(
        _impl->result.qe - _impl->result.qb + 1,
        _impl->ref.get() + _impl->result.qb,
        _impl->result.te - _impl->result.tb + 1,
        _impl->alt.get() + _impl->result.tb,
        5, _impl->mat, _impl->gapo, _impl->gape,
        _impl->altlen,
        &_impl->cigar_len, &_impl->cigar);

    _impl->valid_result = true;
}
