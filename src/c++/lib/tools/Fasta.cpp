
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
#include <cstdio>
#include <limits>
#include <map>
#include <fstream>
#include <mutex>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <boost/algorithm/string.hpp>

extern "C" {

// GCC warns us about some things in htslib here. We don't care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"

#include <htslib/faidx.h>

#pragma GCC diagnostic pop
}

//#define DEBUG_FASTAFILE

class MMappedFastaFile {
public:
    explicit MMappedFastaFile(std::string const & _filename) :
        filename(_filename)
    {
        struct stat st;
        stat(_filename.c_str(), &st);
        filesize = st.st_size;
        fd = open(_filename.c_str(), O_RDONLY, 0);
        assert(fd != -1);
        base = (uint8_t *) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
        assert(base != MAP_FAILED);

        auto fai_fp = fopen((filename + ".fai").c_str(), "r");
        if(!fai_fp)
        {
            assert(fai_build(filename.c_str()) == 0);
        }
        else
        {
            fclose(fai_fp);
        }

        std::ifstream index_in(filename + ".fai");
        while(index_in.good())
        {
            std::string line;
            std::getline(index_in, line);
            std::vector<std::string> parts;
            stringutil::split(line, parts, "\n\t", false);
            if(parts.size() == 5)
            {
                const std::string contig = parts[0];
                index_entry ientry;
                sscanf(parts[1].c_str(), "%zu", &ientry.length);
                sscanf(parts[2].c_str(), "%zu", &ientry.start_offset);
                sscanf(parts[3].c_str(), "%zu", &ientry.chars_per_line);
                sscanf(parts[4].c_str(), "%zu", &ientry.bytes_per_line);
                fai[contig] = ientry;
            }
            else if(parts.size() > 0)
            {
                error("invalid fai line %s in %s", line.c_str(), (filename + ".fai").c_str());
            }
        }
    }

    ~MMappedFastaFile()
    {
        const int rc = munmap(base, filesize);
        assert(rc == 0);
        close(fd);
    }


    std::string get(std::string const & contig, size_t start, size_t len)
    {
        auto ientry = fai.find(contig);
        if(ientry == fai.end())
        {
            error("Contig %s is not known", contig.c_str());
        }

        if(start >= ientry->second.length)
        {
            return "";
            /* TODO: handle downstream */
            /* error("Position %zu is past the end of contig %s", start, contig.c_str()); */
        }

        const size_t line = start / ientry->second.chars_per_line;
        size_t offset_in_line = start % ientry->second.chars_per_line;
        size_t chars_in_line = ientry->second.chars_per_line - offset_in_line;
        size_t line_offset = ientry->second.start_offset + line*ientry->second.bytes_per_line;

        chars_in_line = std::min(chars_in_line, ientry->second.length - start);

        std::string result;
        while(len > 0)
        {
            const size_t to_read = std::min(len, chars_in_line);
            result += std::string((const char *)(base + line_offset + offset_in_line), to_read);
            len -= to_read;
            start += to_read;
            offset_in_line = 0;
            line_offset += ientry->second.bytes_per_line;
            if(start >= ientry->second.length)
            {
                break;
            }
            chars_in_line = std::min(ientry->second.chars_per_line, ientry->second.length - start);
        }
        return result;
    }

private:
    typedef struct _index_entry
    {
        size_t length;
        size_t start_offset;
        size_t chars_per_line;
        size_t bytes_per_line;
    } index_entry;

    std::string filename;
    int fd;
    size_t filesize;
    uint8_t * base;

    std::map<std::string, index_entry> fai;
};


class FastaFileCache
{
public:
    std::shared_ptr<MMappedFastaFile> operator() (std::string const & filename)
    {
        std::lock_guard<std::mutex> write_lock(write_mutex);
        auto it = cache.find(filename);
        if(it == cache.end())
        {
#ifdef DEBUG_FASTAFILE
            std::cerr << "fastafile: adding " << filename << std::endl;
#endif
            it = cache.emplace(filename, std::make_shared<MMappedFastaFile>(filename)).first;

#ifdef DEBUG_FASTAFILE
            std::cerr << "fastafile: added " << filename << std::endl;
        }
        else
        {
            std::cerr << "fastafile: reused " << filename << std::endl;
#endif
        }
        return it->second;
    }
private:
    static std::map<std::string, std::shared_ptr<MMappedFastaFile>> cache;
    static std::mutex write_mutex;
};

std::map<std::string, std::shared_ptr<MMappedFastaFile>> FastaFileCache::cache;
std::mutex FastaFileCache::write_mutex;

struct FastaFileImpl
{
    explicit FastaFileImpl(const char * _filename) :
        filename(_filename)
    {
        FastaFileCache cache;
        contigs = cache(filename);
    }

    ~FastaFileImpl() {}

    std::string filename;
    std::shared_ptr<MMappedFastaFile> contigs;
};

FastaFile::FastaFile() {
    _impl = NULL;
}

FastaFile::FastaFile(const char * filename)
{
    _impl = new FastaFileImpl(filename);
}

FastaFile::~FastaFile()
{
    if(_impl) {
        delete _impl;
    }
}

std::string FastaFile::getFilename() const
{
    if(!_impl) return "";
    return _impl->filename;
}

FastaFile::FastaFile(FastaFile const & rhs)
{
    _impl = new FastaFileImpl(rhs._impl->filename.c_str());
}

FastaFile & FastaFile::operator=(FastaFile const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    if(_impl) delete _impl;
    _impl = new FastaFileImpl(rhs._impl->filename.c_str());
    return *this;
}


std::string FastaFile::query(std::string const & location) const
{
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;

    stringutil::parsePos(location, chr, start, end);
    return query(chr.c_str(), start, end);
}

std::string FastaFile::query(const char * chr, int64_t start, int64_t end) const
{
    if(end < start)
    {
        return "";
    }

    start = std::max(start, 0l);

    const int64_t requested_length = end - start + 1;
    if(requested_length <= 0)
    {
        return "";
    }

    std::string result = _impl->contigs->get(chr, start, requested_length);
    boost::to_upper(result);
    return result;
}
