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
 * \brief Interval buffer implementation
 *
 * \file IntervalBuffer.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "helpers/IntervalBuffer.hh"

#include <algorithm>
#include <map>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "Error.hh"

#include "helpers/IntervalList.hh"

namespace intervals
{

struct IntervalBuffer::IntervalBufferImpl
{
    IntervalBufferImpl() {}
    IntervalBufferImpl(IntervalBufferImpl const & rhs) : lanes(rhs.lanes) {}
    typedef IntervalList<interval> ivmap_t;

    std::vector<ivmap_t> lanes;
};

/** tracks intervals over a number of lanes */
IntervalBuffer::IntervalBuffer()
{
    _impl = new IntervalBufferImpl();
}

IntervalBuffer::~IntervalBuffer()
{
    delete _impl;
}

IntervalBuffer::IntervalBuffer(IntervalBuffer const & rhs)
{
    _impl = new IntervalBufferImpl(*rhs._impl);
}

IntervalBuffer & IntervalBuffer::operator=(IntervalBuffer const & rhs)
{
    if(&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new IntervalBufferImpl(*rhs._impl);
    return *this;
}


/**
 * @brief Add an interval to a lane
 *
 * @param start interval coordinates
 * @param end interval coordinates
 * @param lane lane to add to
 * @return Interval identifier
 */
void IntervalBuffer::addInterval(int64_t start, int64_t end, size_t lane)
{
    if (start > end)
    {
        return;
    }
    if (_impl->lanes.size() <= lane)
    {
        _impl->lanes.resize(lane + 1);
    }
    _impl->lanes[lane].add(interval(start, end));
}

/**
 * @brief Advance buffer, discarding all intervals with start < to
 *
 * @param to interval minimum end position
 */
void IntervalBuffer::advance(int64_t to)
{
    if(to < 0)
    {
        _impl->lanes.clear();
        return;
    }

    for (size_t lane = 0; lane < _impl->lanes.size(); ++lane)
    {
        _impl->lanes[lane].remove_to(to-1);
    }
}

/**
 * @brief Check if interval is fully covered in a given lane
 */
bool IntervalBuffer::isCovered(int64_t start, int64_t end, size_t lane) const
{
    if (lane >= _impl->lanes.size())
    {
        return false;
    }

    // intervals of zero length count as covered
    if (end < start)
    {
        return true;
    }

    std::list<interval> is;
    _impl->lanes[lane].get(start, end, is);
    if (is.size() != 1)
    {
        // if we overlap with more than one interval, then there must be a gap
        return false;
    }

    interval & it = is.front();
    return it.start <= start && it.end >= end ? true : false;
}

/**
 * @brief Check if interval is partially covered in a given lane
 */
bool IntervalBuffer::hasOverlap(int64_t start, int64_t end, size_t lane) const
{
    if (lane >= _impl->lanes.size())
    {
        return false;
    }

    // intervals of zero length count as covered
    if (end < start)
    {
        return true;
    }

    interval it = _impl->lanes[lane].query(start, end);
    return it.start >= 0 && it.end >= 0 && it.end - it.start + 1 > 0;
}

}
