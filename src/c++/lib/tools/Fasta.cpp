
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
 *  \brief Fasta wrapper implementation
 *  
 * \file Fasta.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Fasta.hh"
#include "Error.hh"
#include "helpers/StringUtil.hh"

#include <boost/filesystem.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <fstream>

#include <boost/algorithm/string.hpp>

extern "C" {

// GCC warns us about some things in htslib here. We don't care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"

#include <htslib/faidx.h>

#pragma GCC diagnostic pop
}

struct FastaFileImpl
{
    FastaFileImpl(const char * _filename) : 
        filename(_filename)
    {
        boost::filesystem::path iname(_filename);
        iname += ".fai";
        if(!boost::filesystem::exists(iname))
        {
            if(fai_build(_filename) < 0)
            {
                error("Cannot index %s", _filename);
            }            
        }
        idx = fai_load(_filename);
        if(!idx)
        {
            error("Cannot load index for %s", _filename);
        }

        // read contig lengths from fai since htslib doesn't expose this
        std::ifstream fai(iname.c_str());
        while(fai.good())
        {
            std::string line;
            std::vector<std::string> v;
            std::getline(fai, line);
            stringutil::split(line, v, "\t");
            // fai entries have 5 columns
            if(v.size() == 5)
            {
                contig_lengths[v[0]] = std::stol(v[1]);
                // std::cerr << v[0] << ":" << std::stol(v[1]) << "\n";
            }
        }
    }

    ~FastaFileImpl()
    {
        fai_destroy(idx);
    }

    faidx_t * idx;
    std::string filename;
    std::map<std::string, int64_t> contig_lengths;
};

FastaFile::FastaFile(const char * filename)
{
    _impl = new FastaFileImpl(filename);
}

FastaFile::~FastaFile()
{
    delete _impl;
}

std::string FastaFile::getFilename() const
{
    return _impl->filename;
}

FastaFile::FastaFile(FastaFile const & rhs)
{
    _impl = new FastaFileImpl(rhs._impl->filename.c_str());
}

FastaFile const & FastaFile::operator=(FastaFile const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new FastaFileImpl(rhs._impl->filename.c_str());
    return *this;
}


std::string FastaFile::query(std::string const & location)
{
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;

    stringutil::parsePos(location, chr, start, end);
    return query(chr.c_str(), start, end);
}

std::string FastaFile::query(const char * chr, int64_t start, int64_t end)
{
    int64_t requested_length = end-start+1;

    auto clen = _impl->contig_lengths.find(chr);
    if(clen != _impl->contig_lengths.end())
    {
        if(start+1 > clen->second)
        {
            return "";
        }
    }
    else
    {
        std::cerr << "[W] FastaFile::query:  Unknown contig length for " << chr << " --  we might read over the end." << "\n";
    }

    if(end < 0 || start < 0)
    {
        requested_length = std::numeric_limits<int64_t>::max();
    }

    if(requested_length <= 0)
    {
        return "";
    }
    int len;
    // faidx_fetch_seq (..., start, end) gets [start, end]
    char * data = faidx_fetch_seq(_impl->idx, chr, start, end, &len);

    if(len <= 0)
    {
        return "";
    }

    std::string result(data, std::min(requested_length, (int64_t)len));
    free(data);
    boost::to_upper(result);
    return result;
}
