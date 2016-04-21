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
 *  \brief String helper functions
 *
 * \file StringUtil.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <string>
#include <vector>
#include <sstream>

namespace stringutil
{
    /**
     * @brief Split a string along separators
     */
    static inline void split(std::string str,
               std::vector<std::string > & result,
               const std::string & seps = " ,",
               bool return_empty=false
              )
    {
        bool has_last = false;
        while (str.size() > 0)
        {
            size_t pos = str.find_first_of(seps);
            if (pos != std::string::npos)
            {
                if(return_empty || pos > 0)
                {
                    result.push_back(str.substr(0, pos));
                }
                str = str.substr(pos + 1);
                has_last = true;
            }
            else
            {
                result.push_back(str);
                str = "";
                has_last = false;
            }
        }
        if (has_last && return_empty)
        {
            result.push_back("");
        }
    }

    /**
     * @brief Join a list of strings / values with a separator
     *
     */
    template <class _strcontainer>
    static inline std::string join(_strcontainer const & c, const std::string & sep)
    {
        std::string result;
        bool first = true;
        for(auto const & x : c)
        {
            if(!first)
            {
                result += sep;
            }
            else
            {
                first = false;
            }
            result += x;
        }
        return result;
    }

    /**
     * @brief Test if string has given suffix
     *
     * @return true if str ends with suffix
     */
    static inline bool endsWith(std::string const & str, std::string const & suffix)
    {
        if(suffix.size() > str.size())
        {
            return false;
        }
        return str.substr(str.size() - suffix.size()) == suffix;
    }

    /**
     * @brief Replace all instances of find with replace in str
     */
    static inline std::string replaceAll(std::string str, std::string const & find, std::string const & replace)
    {
        size_t start_pos = 0;
        while((start_pos = str.find(find, start_pos)) != std::string::npos)
        {
            str.replace(start_pos, find.length(), replace);
            start_pos += replace.length();
        }
        return str;
    }

    /**
     * @brief Replace all instances of find with replace in str
     */
    static inline void replaceAllInplace(std::string & str, std::string const & find, std::string const & replace)
    {
        size_t start_pos = 0;
        while((start_pos = str.find(find, start_pos)) != std::string::npos)
        {
            str.replace(start_pos, find.length(), replace);
            start_pos += replace.length();
        }
    }

    /**
     * Format genomic coordinates
     */
    static inline std::string formatPos(const char * chr, int64_t pos = -1, int64_t end = -1)
    {
        std::stringstream ss;

        if(chr)
        {
            ss << chr;
        }
        if (pos >= 0)
        {
            ss << ":";
            ss << (pos+1);
            if (end >= 0)
            {
                ss << "-" << (end+1);
            }
        }

        return ss.str();
    }

    /**
     * Format genomic coordinates
     */
    static inline std::string formatPos(std::string const & chr, int64_t pos = -1, int64_t end = -1)
    {
        return formatPos(chr.c_str(), pos, end);
    }


    /**
     * @brief Parse coordinates
     *
     * @param input input string, e.g. "chr1:1,000-2000"
     * @param chr the contig name, e.g. "chr1"
     * @param start the start, e.g. 999
     * @param end the end, e.g. 1999
     *
     * Returned coordinates are 0-based, input coordinates are 1-based.
     *
     */
    static inline void parsePos(std::string input, std::string & chr, int64_t & start, int64_t & end)
    {
        std::vector<std::string> spl;
        split(input, spl, " :-");

        if(spl.size() >= 1)
        {
            chr = spl[0];
        }
        if(spl.size() >= 2)
        {
            start = std::stoll(replaceAll(spl[1], ",", "")) - 1;
        }
        if(spl.size() >= 3)
        {
            end = std::stoll(replaceAll(spl[2], ",", "")) - 1;
        }
    }
}
