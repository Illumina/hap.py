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
 *  \brief Wrapper for htslib faidx
 *
 * \file Fasta.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <list>
#include <string>

struct FastaFileImpl;
class FastaFile
{
public:
    FastaFile();
	explicit FastaFile(const char * filename);

	~FastaFile();
	FastaFile(FastaFile const &);
	FastaFile & operator=(FastaFile const &);

	std::string getFilename() const;

	std::string query(std::string const & location) const;
	std::string query(const char * chr, int64_t start, int64_t end) const;

	/**
	 * return the size of a contig. Sum of all contig sizes if contig == 0
	 * @param contig name of contig, or empty
	 * @return number of bases in the contig
	 */
	size_t contigSize(std::string const & contig = "") const;

	/**
	 * return the non-N padded size of a contig. This is calculated
	 * as the size of the contig minus any N's at the beginning or
	 * at the end.
	 *
	 * Returns the sum of all contig sizes if contig == 0
	 * @param contig name of contig, or empty
	 * @return number of bases in the contig
	 */
	size_t contigNonNSize(std::string const & contig = "") const;

	/**
	 * @return all contig names
	 */
	std::list<std::string> getContigNames() const;
private:
	FastaFileImpl * _impl;
};
