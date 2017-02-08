
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
#include <errno.h>
#include <string.h>

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
        filesize = (size_t) st.st_size;
        fd = open(_filename.c_str(), O_RDONLY, 0);
        assert(fd != -1);
#if __APPLE__
        base = (uint8_t *) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, fd, 0);
#else
        base = (uint8_t *) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
#endif
        if(base == MAP_FAILED)
        {
            const int err = errno;
            error("Cannot mmap %s (errno=%i / %s) -- do you have enough memory available?",
                  _filename.c_str(), err, strerror(err));
        }

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

                // determine non-N length; we don't do this because it's pretty slow.
//                size_t offset = ientry.start_offset;
//                const size_t line_number = (ientry.length - 1) / ientry.chars_per_line;
//                size_t offset_in_line = (ientry.length - 1) % ientry.chars_per_line;
//                size_t offset_end = ientry.start_offset + line_number*ientry.bytes_per_line + offset_in_line + 1;
//                size_t ns = 0;
//
//                while(offset < offset_end)
//                {
//                    if(std::tolower(base[offset - 1]) == 'n')
//                    {
//                        ++ns;
//                    }
//                    ++offset;
//                }
//                fai[contig].non_n_length = fai[contig].length - ns;

                // length trimming off N's at start and end
                size_t pos = 0;
                size_t ns_at_start = 0;
                bool done = false;

                // start of contig
                while(pos < ientry.length && !done)
                {
                    const std::string s = get(contig, pos, 10000);

                    for(size_t j = 0; j < s.size(); ++j)
                    {
                        if(std::tolower(s.at(j)) == 'n')
                        {
                            ++ns_at_start;
                        }
                        else
                        {
                            done = true;
                            break;
                        }
                    }

                    pos += s.size();
                }

                size_t ns_at_end = 0;
                // check if we had all Ns
                if(ns_at_start < ientry.length)
                {
                    const size_t line_number = (ientry.length - 1) / ientry.chars_per_line;
                    size_t offset_in_line = (ientry.length - 1) % ientry.chars_per_line;
                    size_t offset = ientry.start_offset + line_number*ientry.bytes_per_line + offset_in_line + 1;

                    // end of contig
                    while(offset > ientry.start_offset)
                    {
                        if (base[offset - 1] == '\n' || base[offset - 1] == '\r')
                        {
                            // skip newlines
                        }
                        else if(std::tolower(base[offset - 1]) == 'n')
                        {
                            ++ns_at_end;
                        }
                        else
                        {
                            break;
                        }
                        --offset;
                    }
                }
                fai[contig].non_n_length = fai[contig].length - (ns_at_start + ns_at_end);
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


    /**
     * return the size of a contig. Sum of all contig sizes if contig == 0
     * @param contig name of contig, or empty
     * @return number of bases in the contig
     */
    size_t contigSize(std::string const & contig) const
    {
        size_t result = 0;
        if(contig.empty())
        {
            for(const auto & ix : fai)
            {
                result += ix.second.length;
            }
        }
        else
        {
            auto k = fai.find(contig);
            if(k != fai.end())
            {
                result += k->second.length;
            }
        }

        return result;
    }

    /**
     * return the non-N padded size of a contig. This is calculated
     * as the size of the contig minus any N's at the beginning or
     * at the end.
     *
     * Returns the sum of all contig sizes if contig == 0
     * @param contig name of contig, or empty
     * @return number of bases in the contig
     */
    size_t contigNonNSize(std::string const & contig) const
    {
        size_t result = 0;
        if(contig.empty())
        {
            for(const auto & ix : fai)
            {
                result += ix.second.non_n_length;
            }
        }
        else
        {
            auto k = fai.find(contig);
            if(k != fai.end())
            {
                result += k->second.non_n_length;
            }
        }

        return result;
    }

    std::list<std::string> getContigNames() const
    {
        std::list<std::string> names;
        for(auto const & c : fai)
        {
            names.push_back(c.first);
        }
        return names;
    }

private:
    typedef struct _index_entry
    {
        size_t length;
        size_t start_offset;
        size_t chars_per_line;
        size_t bytes_per_line;
        size_t non_n_length;
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

    start = std::max(start, (int64_t )0l);

    const int64_t requested_length = end - start + 1;
    if(requested_length <= 0)
    {
        return "";
    }

    std::string result = _impl->contigs->get(chr, (size_t) start, (size_t) requested_length);
    boost::to_upper(result);
    return result;
}

/**
 * return the size of a contig. Sum of all contig sizes if contig == 0
 * @param contig name of contig, or empty
 * @return number of bases in the contig
 */
size_t FastaFile::contigSize(std::string const & contig) const
{
    return _impl->contigs->contigSize(contig);
}

/**
 * return the non-N padded size of a contig. This is calculated
 * as the size of the contig minus any N's at the beginning or
 * at the end.
 *
 * Returns the sum of all contig sizes if contig == 0
 * @param contig name of contig, or empty
 * @return number of bases in the contig
 */
size_t FastaFile::contigNonNSize(std::string const & contig) const
{
    return _impl->contigs->contigNonNSize(contig);
}

/**
 * @return all contig names
 */
std::list<std::string> FastaFile::getContigNames() const
{
    return _impl->contigs->getContigNames();
}
