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
 * \brief Buffer for sorted intervals
 *
 *
 * \file IntervalBuffer.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <list>
#include <utility>

namespace intervals
{

class IntervalBuffer
{
public:
    /** tracks intervals over a number of lanes */
    IntervalBuffer();
    IntervalBuffer(IntervalBuffer const & rhs);
    ~IntervalBuffer();
    IntervalBuffer & operator=(IntervalBuffer const & rhs);

    /**
     * @brief Add an interval to a lane
     *
     * @param start interval coordinates
     * @param end interval coordinates
     * @param lane lane to add to
     */
    void addInterval(int64_t start, int64_t end, size_t lane);

    /**
     * @brief Advance buffer, discarding all intervals with end < to
     *
     * @param to interval minimum end position; pass -1 to clear buffer
     */
    void advance(int64_t to);

    /**
     * @brief Check if interval is fully covered in a given lane
     */
    bool isCovered(int64_t start, int64_t end, size_t lane) const;

    /**
     * @brief Check if interval is partially covered in a given lane
     */
    bool hasOverlap(int64_t start, int64_t end, size_t lane) const;

    /**
     * Intersect two lanes and return the size of the intersection
     * @return the number of positions shared by l1 and l2
     */
    size_t intersectLanes(const size_t l1, const size_t l2) const;

private:
    struct IntervalBufferImpl;
    IntervalBufferImpl * _impl;
};

}
