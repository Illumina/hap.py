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
 * Some bit helpers
 *
 * \file Popcount.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cstdlib>

#ifdef __GNUC__

#define popcount64 __builtin_popcountl

// http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
static inline uint64_t next_permutation(uint64_t v)
{
    uint64_t w; // next permutation of bits
    uint64_t t = v | (v - 1); // t gets v's least significant 0 bits set to 1

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzl(v) + 1));
    return w;
}

#else

// https://en.wikipedia.org/wiki/Hamming_weight
static inline int popcount64(uint64_t x)
{
    const uint64_t m1  = 0x5555555555555555; //binary: 0101...
    const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
    const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

// http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
static inline uint64_t next_permutation(uint64_t v)
{
    uint64_t t = (v | (v - 1)) + 1;
    w = t | ((((t & -t) / (v & -v)) >> 1) - 1);
    return w;
}

#endif

static inline uint64_t rand64()
{
    // assumes RAND_MAX < 0x10000
    return ( rand() & 0xffff ) |
           ( (uint64_t )( rand() & 0xffff ) << 16) |
           ( (uint64_t )( rand() & 0xffff ) << 32) |
           ( (uint64_t )( rand() & 0xffff ) << 48);
}
