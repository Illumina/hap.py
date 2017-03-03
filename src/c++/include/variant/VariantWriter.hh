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
 * Simple Variant-to-file writer
 *
 * \file VariantWriter.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Variant.hh"

namespace variant
{

/**
 * @brief Write variants to file
 */
struct VariantWriterImpl;
class VariantWriter {
public:
    VariantWriter(const char * filename, const char * reference);

    VariantWriter(VariantWriter const &);
    ~VariantWriter();
    VariantWriter & operator=(VariantWriter const &);

    /** write format fields apart from GT */
    void setWriteFormats(bool write_fmts=false);
    bool getWriteFormats() const;

    /**
     * @brief Get header from VariantReader
     *
     * Optionally add only contig names, but not format and ID fields
     */
    void addHeader(VariantReader &, bool drop_formats_info=false);

    /**
     * @brief Add line to VCF header
     *
     * @param headerline the line to add
     */
    void addHeader(const char * headerline);

    /**
     * @brief Add a sample to the header
     * @param sname  name of sample to write
     * @return sample number to retrieve records in get()
     */
    int addSample(const char * sname);

    /**
     * @brief write a variant to a file
     *
     * @param var the variant records to write
     */
    void put(Variants const & var);

private:
    VariantWriterImpl * _impl;
};

}
