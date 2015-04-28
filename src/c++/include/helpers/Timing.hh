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
 * \brief Helper to get CPU time
 *
 * See http://stackoverflow.com/questions/275004/c-timer-function-to-provide-time-in-nano-seconds 
 *
 * \file Timing.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <chrono>

struct CPUClock
{
    typedef unsigned long long                 rep;
    typedef std::ratio<1, 2000000000>          period; // My machine is 2 GHz
    typedef std::chrono::duration<rep, period> duration;
    typedef std::chrono::time_point<CPUClock>     time_point;
    static const bool is_steady =              true;

    static time_point now() noexcept
    {
        unsigned lo, hi;
        asm volatile("rdtsc" : "=a" (lo), "=d" (hi));
        return time_point(duration(static_cast<rep>(hi) << 32 | lo));
    }
};


#define TIMEIT(result, N, ts, ...) do {                                             \
    typedef std::chrono::duration<double, typename CPUClock::period> Cycle;         \
    using std::chrono::duration_cast;                                               \
                                                                                    \
    auto t0 = CPUClock::now();                                                      \
    for(int n = 0; n < N; ++n)                                                      \
    {                                                                               \
        ts timeme(__VA_ARGS__);                                                     \
    }                                                                               \
    auto t1 = CPUClock::now();                                                      \
    auto t2 = CPUClock::now();                                                      \
    for(int n = 0; n < N; ++n)                                                      \
    {                                                                               \
        ts timeme(__VA_ARGS__);                                                     \
        timeme();                                                                   \
    }                                                                               \
    auto t3 = CPUClock::now();                                                      \
    auto ticks_per_iter = Cycle((t3-t2) - (t1-t0))/N;                               \
    result = ticks_per_iter.count();                                                \
    } while(0)
