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
 * \file Klib.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Alignment.hh"

struct KlibAlignmentImpl;
class KlibAlignment : public Alignment
{
public:
    KlibAlignment();
    ~KlibAlignment();

    void setParameters(AlignmentParameters const & ap);
    void getParameters(AlignmentParameters & ap);
    
    /**
     * @brief set target sequence
     */
    void setRef(const char * seq);
    /*
     * @brief set query sequence
     */
    void setQuery(const char * seq);
    /**
     * @brief Get the alignment score
     */
    int getScore();
    /**
     * @brief Get a cigar string
     */
    void getCigar(int &r0, int & r1, int &a0, int & a1, std::string &cig);
    void getCigar(int &r0, int & r1, int &a0, int & a1, int & n_cigar, uint32_t *& cigar);

protected:
    virtual void update();
	KlibAlignmentImpl * _impl;
};

