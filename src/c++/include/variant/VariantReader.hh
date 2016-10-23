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
 * Reader for getting Variants from File
 *
 * \file VariantReader.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"

#include <list>
#include <string>

namespace variant
{
/**
 * @brief read variants from a file
 * @details Wrapper for synced_bcf_reader
 */
struct VariantReaderImpl;
class VariantReader {
public:
    VariantReader();

    VariantReader(VariantReader const &);
    ~VariantReader();
    VariantReader const & operator=(VariantReader const &);

    /**
     * @brief Apply filters in the input VCF(s)
     *
     * Filters can be enabled / disabled globally (sample == -1)
     * or for each sample individually.
     *
     */
    void setApplyFilters(bool filters=false, int sample = -1);
    bool getApplyFilters(int sample = -1) const;

    /**
     * @brief Return homref/no-calls
     *
     */
    void setReturnHomref(bool homref=true);
    bool getReturnHomref() const;

    /**
     * @brief change GTs on chrX/Y to be diploid for matching
     *
     */
    void setFixChrXGTs(bool fix=true);

    /**
     * @brief Validate reference alleles
     *
     */
    void setValidateRef(const char * ref_fasta, bool validate=true);
    bool getValidateRef() const;

    /**
     * @brief Interface to htslib regions functionality
     * @param regions regions string, see synced_bcf_reader.h
     * @param isFile True if regions is a file
     *
     * Must be called before addSample!
     *
     */
    void setRegions(const char * regions, bool isFile);

    /**
     * @brief Interface to htslib targets functionality
     * @param targets targets string, see synced_bcf_reader.h
     * @param isFile True if targets is a file
     *
     * Must be called before addSample!
     *
     */
    void setTargets(const char * targets, bool isFile);

    /**
     * @brief Add a sample to read from
     * @param filename  file name
     * @param sname  name of sample to read from (or "" to use first sample in file / use "*" to add all samples in a file)
     * @return sample number to retrieve records in the calls vector in current()
     */
    int addSample(const char * filename, const char * sname);

    /** return a list of filename/sample pairs */
    void getSampleList(std::list< std::pair<std::string, std::string> > & files);

    /**
     * @brief Rewind / set region to read
     *
     * @param chr chromosome/contig name
     * @param startpos start position on chr (or -1 for entire chr)
     */
    void rewind(const char * chr=NULL, int64_t startpos=-1);

    /**
     * @brief Return next variant and advance
     *
     * @param v Variant record to populate
     */
    Variants & current();

    /**
     * @brief Advance one position
     * @return true if a variant was retrieved, false otherwise
     */
    bool advance();

    /**
     * @brief Insert a Variants record into the stream
     * @details The next record returned will be this one
     *
     * @param s Variants to put on top of stack
     * @param back enqueue at back or front of buffer
     */
    void enqueue(Variants const & v, bool back=false);

private:
    friend class VariantWriter;
    VariantReaderImpl * _impl;
};

}
