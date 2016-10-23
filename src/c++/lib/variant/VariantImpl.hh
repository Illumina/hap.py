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
 * \brief Implementation header for VariantReader/Writer
 *
 * \file VariantImpl.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Variant.hh"
#include "Fasta.hh"
#include "Error.hh"

#include <cassert>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <vector>
#include <limits>
#include <cmath>

#include "helpers/BCFHelpers.hh"
#include "helpers/StringUtil.hh"

#include <boost/algorithm/string.hpp>


namespace variant
{

typedef struct _SampleInfo
{
    std::string filename;
    std::string sample;
    int ireader;
    int isample;
    int end_in_vcf;
    std::vector<int> allele_map;
} SampleInfo;


/** Private data for VariantReader */
struct VariantReaderImpl
{
    VariantReaderImpl()
    {
        files = bcf_sr_init();
        files->require_index = 1;
        files->collapse = COLLAPSE_NONE;
        regions = "";
        targets = "";
        regionsFile = false;
        targetsFile = false;
        applyFilters = false;
        returnHomref = true;
        validateRef = false;
        fix_chrX = false;
    }

    ~VariantReaderImpl()
    {
        bcf_sr_destroy(files);
    }

    bcf_srs_t *files;

    // keep filenames for each sample
    std::vector<SampleInfo> samples;

    // mapping of file names to readers
    std::map<std::string, int> filename_mapping;

    // setRegions / setTargets
    std::string regions;
    bool regionsFile;
    std::string targets;
    bool targetsFile;

    // apply filters
    bool applyFilters;
    std::vector<bool> applyFiltersPerSample;

    // return homref / no-call variants
    bool returnHomref;

    // validate ref alleles
    FastaFile ref;
    bool validateRef;

    // we buffer variant output
    std::list<Variants> buffered_variants;

    bool fix_chrX;
};

struct VariantWriterImpl
{
    VariantWriterImpl(const char * fname, const char * refname) :
        write_formats(false), filename(fname), referencename(refname),
        header_done(false), reference(refname)
    {
        const char * mode = "wu";

        if(stringutil::endsWith(fname, ".vcf.gz"))
        {
            mode = "wz";
        }
        else if(stringutil::endsWith(fname, ".bcf"))
        {
            mode = "wb";
        }

        if(strlen(fname) > 0 && fname[0] == '-')
        {
            fp = hts_open("-", mode);
        }
        else
        {
            fp = hts_open(fname, mode);
        }

        hdr = bcf_hdr_init("w");
        bcfhelpers::bcfHeaderHG19(hdr);
        rec = bcf_init1();
    }

    ~VariantWriterImpl()
    {
        bcf_destroy1(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
    }

    void writeHeader();

    std::vector< std::string > header_lines;

    bool write_formats;

    std::string filename;
    std::string referencename;
    std::vector<std::string> samples;

    // internal
    bool header_done;

    htsFile * fp;
    bcf_hdr_t *hdr;
    bcf1_t *rec;

    FastaFile reference;
};

} // namespace variant
