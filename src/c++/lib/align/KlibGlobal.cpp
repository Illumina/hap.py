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
 *  \brief Implementation: global alignment via klib
 *
 *
 * \file KlibGlobal.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "KlibGlobal.hh"

#include "KlibImpl.hh"
#include "Error.hh"

// see ksw.c
#define MINUS_INF -0x40000000

void KlibGlobalAlignment::update()
{
	_impl->result.qb = 0;
	_impl->result.qe = _impl->reflen-1;
	_impl->result.tb = 0;
	_impl->result.te = _impl->altlen-1;

	if(_impl->cigar)
	{
		free(_impl->cigar);
		_impl->cigar = NULL;
		_impl->cigar_len = 0;
	}

	_impl->result.score = ksw_global(
		_impl->result.qe - _impl->result.qb + 1,
		_impl->ref.get() + _impl->result.qb,
		_impl->result.te - _impl->result.tb + 1,
		_impl->alt.get() + _impl->result.tb,
		5, _impl->mat, _impl->gapo, _impl->gape,
		std::max(_impl->reflen, _impl->altlen),
		&_impl->cigar_len, &_impl->cigar);

	if(_impl->result.score <= MINUS_INF)
	{
		error("Failed to globally align.");
	}

	_impl->valid_result = true;
}
