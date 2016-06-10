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
 * \brief Location info data structure
 *
 * \details Store values associated with integer locations
 *
 *
 * \file LocationInfo.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <list>
#include <utility>
#include <type_traits>

#include "helpers/IntervalList.hh"

namespace intervals
{

/** work only with types we can memcpy */
template<class value_t,
         class T = typename std::enable_if< std::is_fundamental<value_t>::value>::type>
struct value_interval : public interval
{
    value_interval() : interval(), values(NULL) {}
    value_interval(int64_t _start, int64_t _end) : interval(_start, _end)
    {
        if (end >= start)
        {
            values = new value_t[end - start + 1];
            memset(values, 0, sizeof(value_t)*(end - start + 1));
        }
        else
        {
            values = NULL;
        }
    }

    value_interval(int64_t _start, int64_t _end, value_t const & val) : interval(_start, _end)
    {
        if (end >= start)
        {
            values = new value_t[end - start + 1];
            for (int64_t i = 0; i < end - start + 1; ++i)
            {
                values[i] = val;
            }
        }
        else
        {
            values = NULL;
        }
    }

    value_interval(value_interval const & rhs) : interval(rhs)
    {
        if (end >= start)
        {
            values = new value_t[end - start + 1];
            memcpy(values, rhs.values, sizeof(value_t)*(end - start + 1));
        }
        else
        {
            values = NULL;
        }
    }
    ~value_interval()
    {
        if(values)
        {
            delete [] values;
        }
    }

    value_interval & operator=(value_interval const & rhs)
    {
        if (&rhs == this)
        {
            return *this;
        }
        interval::operator=(rhs);
        if (values)
        {
            delete [] values;
        }
        if (end >= start)
        {
            values = new value_t[end - start + 1];
            memcpy(values, rhs.values, sizeof(value_t)*(end - start + 1));
        }
        else
        {
            values = NULL;
        }
        return *this;
    }


    // merge two intervals
    void merge(interval const & rhs)
    {
        int64_t new_start = std::min(rhs.start, start);
        int64_t new_end = std::max(rhs.end, end);
        value_t * new_values = NULL;

        if(new_end >= new_start)
        {
            new_values = new value_t[new_end - new_start + 1];
            memset(new_values, 0, sizeof(value_t)*(new_end - new_start + 1));

            if (values)
            {
                memcpy(new_values + start - new_start, values, sizeof(value_t)*(end - start + 1));
            }
            if (static_cast<value_interval const &>(rhs).values)
            {
                memcpy(new_values + rhs.start - new_start,
                       static_cast<value_interval const &>(rhs).values, sizeof(value_t)*(rhs.end - rhs.start + 1));
            }
        }
        if(values)
        {
            delete [] values;
        }
        values = new_values;
        start = new_start;
        end = new_end;
    }

    void resize(int64_t new_start, int64_t new_end)
    {
        value_t * new_values = NULL;

        if(new_end >= new_start)
        {
            new_values = new value_t[new_end - new_start + 1];
            memset(new_values, 0, sizeof(value_t)*(new_end - new_start + 1));

            if (values)
            {
                int64_t copy_from = std::max(start, new_start);
                int64_t copy_len = std::min(new_end, end) - copy_from + 1;
                if (copy_len > 0)
                {
                    memcpy(new_values + copy_from - new_start,
                           values + copy_from - start,
                           sizeof(value_t)*copy_len);
                }
            }
        }
        if(values)
        {
            delete [] values;
        }
        values = new_values;
        start = new_start;
        end = new_end;
    }

    value_t * values;
};

/** Store values associated with positive integer locations
 *  Storage assumes mostly dense values (i.e we will
 *  have a value for most locations)
 */
template<class value_t>
class LocationInfo
{
public:
    LocationInfo() {}

    /**
     * Set values in a range
     *
     * @param val the value
     * @param start the start position
     * @param end the end position, or -1 to use the smae as start
     */
    void set(value_t const & val, int64_t start, int64_t end=-1)
    {
        intervals.add(interval_t(start, end, val));
    }

    /**
     * Reset / remove values up to and including position end
     *
     * @param end the end position, or -1 to use the smae as start
     */
    void reset_to(int64_t end = -1)
    {
        intervals.remove_to(end);
    }

    /**
     * Reset / remove values from and including position start
     *
     * @param start the start position, or -1 for all values
     */
    void reset_from(int64_t start = -1)
    {
        intervals.remove_from(start);
    }

    /** query sum over a range */
    size_t querySum(int64_t start, int64_t end, value_t & output)
    {
        output = value_t(0);
        size_t count = 0;
        std::list<interval_t> ivl;
        intervals.get(start, end, ivl);

        auto it = ivl.begin();
        while(start <= end)
        {
            while(it != ivl.end() && it->end < start)
            {
                ++it;
            }
            if(it == ivl.end())
            {
                break;
            }
            start = std::max(start, it->start);
            output += (static_cast<interval_t const &>(*it)).values[start - it->start];
            ++start;
            ++count;
        }
        return count;
    }

    /** query mean value over a range */
    size_t queryMean(int64_t start, int64_t end, value_t & output)
    {
        size_t n = querySum(start, end, output);
        if (n > 0)
        {
            output /= n;
        }
        return n;
    }

    void dump(std::ostream & o)
    {
        for(auto & x : intervals.getIntervals())
        {
            o << x.start << "-" << x.end << "; ";
        }
    }

private:
    typedef value_interval < value_t > interval_t;
    IntervalList<interval_t> intervals;
};

}

