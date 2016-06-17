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
 * \brief VCF checker
 *
 * \details check if vcf can be happily turned into a bcf
 *
 * \file validatevcf.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"

#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>
#include <memory>

#include <htslib/synced_bcf_reader.h>
#include <helpers/BCFHelpers.hh>

#include "helpers/Roc.hh"

// error needs to come after boost headers.
#include "Error.hh"

/* #define DEBUG_VCFCHECK */

namespace WARNING {
    enum WARNING  {
        REFPADDING,
        OVERLAP,
        SIZE
    };
}

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;

    std::string file;

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;
    int64_t rlimit = -1;

    int64_t message = -1;

    bool apply_filters = false;
    bool strict_homref = false;
    bool all_warnings = false;

    try
    {
        try
        {
            // Declare the supported options.
            po::options_description desc("Allowed options");
            desc.add_options()
                ("help,h", "produce help message")
                ("version", "Show version")
                ("input-file", po::value<std::string>(), "The input file")
                ("location,l", po::value<std::string>(), "Start location.")
                ("limit-records", po::value<int64_t>(), "Maximum number of records to process")
                ("message-every", po::value<int64_t>(), "Print a message every N records.")
                ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
                ("strict-homref,H", po::value<bool>(), "Be strict about hom-ref assertions (i.e. don't allow these to overlap).")
                ("all-warnings,W", po::value<bool>(), "Show all warnings, not just the first instance.")
            ;

            po::positional_options_description popts;
            popts.add("input-file", -1);

            po::options_description cmdline_options;
            cmdline_options
                .add(desc)
            ;

            po::variables_map vm;

            po::store(po::command_line_parser(argc, argv).
                      options(cmdline_options).positional(popts).run(), vm);
            po::notify(vm);

            if (vm.count("version"))
            {
                std::cout << "vcfcheck version " << HAPLOTYPES_VERSION << "\n";
                return 0;
            }

            if (vm.count("help"))
            {
                std::cout << desc << "\n";
                return 1;
            }

            if (vm.count("input-file"))
            {
                file = vm["input-file"].as< std::string >();
            }

            if (vm.count("location"))
            {
                stringutil::parsePos(vm["location"].as< std::string >(), chr, start, end);
            }

            if (vm.count("limit-records"))
            {
                rlimit = vm["limit-records"].as< int64_t >();
            }

            if (vm.count("message-every"))
            {
                message = vm["message-every"].as< int64_t >();
            }

            if (vm.count("apply-filters"))
            {
                apply_filters = vm["apply-filters"].as< bool >();
            }

            if (vm.count("strict-homref"))
            {
                strict_homref = vm["strict-homref"].as< bool >();
            }

            if (vm.count("all-warnings"))
            {
                all_warnings = vm["all-warnings"].as< bool >();
            }

            if(file.size() == 0)
            {
                std::cerr << "Please specify one input file / sample.\n";
                return 1;
            }
        }
        catch (po::error & e)
        {
            std::cerr << e.what() << "\n";
            return 1;
        }

        bcf_srs_t * reader = bcf_sr_init();
        reader->collapse = COLLAPSE_NONE;
        if(chr.empty())
        {
            reader->require_index = 0;
            reader->streaming = 1;
        }
        else
        {
            reader->require_index = 1;
            reader->streaming = 0;
        }

        if (!bcf_sr_add_reader(reader, file.c_str()))
        {
            error("Failed to open or file not indexed: %s\n", file.c_str());
        }

        if(!chr.empty())
        {
            int success = 0;
            if(start < 0)
            {
                success = bcf_sr_seek(reader, chr.c_str(), 0);
            }
            else
            {
                success = bcf_sr_seek(reader, chr.c_str(), start);
                std::cerr << "starting at " << chr << ":" << start << "\n";
            }
            if(success < 0)
            {
                error("Cannot seek to %s:%i", chr.c_str(), start);
            }
        }

        bcf_hdr_t * hdr = reader->readers[0].header;

        if(bcf_hdr_nsamples(hdr) < 1)
        {
            error("Input file has no samples. Hap.py will not like that.");
        }

        int nl = 1;
        int64_t rcount = 0;
        std::string current_chr;

        int64_t previous_end = -1;
        // this could be improved using and interval list
        // currently, we only check directly adjacent overlaps
        std::unique_ptr<int[]> previous_allele_count(new int[bcf_hdr_nsamples(hdr)]);
        memset(previous_allele_count.get(), 0, sizeof(int)*bcf_hdr_nsamples(hdr));

        int has_warned[WARNING::SIZE];
        memset(has_warned, 0, sizeof(int)*WARNING::SIZE);

        int n_nonref = 0;

        while(nl)
        {
            nl = bcf_sr_next_line(reader);
            if (nl <= 0)
            {
                break;
            }
            if(!bcf_sr_has_line(reader, 0))
            {
                continue;
            }
            bcf1_t *line = reader->readers[0].buffer[0];

            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }
            const std::string vchr = bcfhelpers::getChrom(hdr, line);
            if(end != -1 && ((!current_chr.empty() && vchr != current_chr) || line->pos > end))
            {
                break;
            }

            if(!current_chr.empty() && vchr != current_chr)
            {
                // chromosome switch
                previous_end = -1;
                memset(previous_allele_count.get(), 0, sizeof(int)*bcf_hdr_nsamples(hdr));
            }

            if(apply_filters)
            {
                bcf_unpack(line, BCF_UN_FLT);

                bool fail = false;
                for(int j = 0; j < line->d.n_flt; ++j)
                {
                    std::string filter = "PASS";
                    int k = line->d.flt[j];

                    if(k >= 0)
                    {
                        filter = bcf_hdr_int2id(hdr, BCF_DT_ID, line->d.flt[j]);
                    }
                    if(filter != "PASS")
                    {
                        fail = true;
                        break;
                    }
                }

                // skip failing
                if(fail)
                {
                    continue;
                }
            }

            bcf_unpack(line, BCF_UN_ALL);

            int64_t vstart = -1, vend= -1;
            bcfhelpers::getLocation(hdr, line, vstart, vend);
            // check
            const int refpadding = bcfhelpers::isRefPadded(line);
            if(refpadding)
            {
#ifdef DEBUG_VCFCHECK
                std::string alleles;
                for(int j = 0; j < line->n_allele; ++j)
                {
                    if(!alleles.empty())
                    {
                        alleles += ",";
                    }
                    alleles += line->d.allele[j];
                }
                std::cerr << "at " << vchr << ":" << vstart << " " << alleles << " refpadding = " << refpadding << "\n";
#endif
                ++vstart;
                if(refpadding > 1)
                {
                    if(all_warnings || !has_warned[WARNING::REFPADDING])
                    {
                        std::cerr << "[W] variant at " << vchr << ":" << vstart << " has more than one base of reference padding \n";
                    }
                    has_warned[WARNING::REFPADDING]++;
                }
            }

            if(line->n_sample != bcf_hdr_nsamples(hdr))
            {
                error("Number of samples in line and header disagrees at %s:%i", vchr.c_str(), vstart);
            }

            bool any_alts = false;
            for(int isample = 0; isample < line->n_sample; ++isample)
            {
                int gt[MAX_GT];
                int ngt;
                bool phased;
                bcfhelpers::getGT(hdr, line, isample, gt, ngt, phased);
                int allele_count = 0;
                for(int g = 0; g < ngt; ++g)
                {
                    if(gt[g] > 0)
                    {
                        const char * alt = line->d.allele[gt[g]];
                        if(*alt == 0 || *alt == '*' || *alt == '<' || *alt == '.')
                        {
                            // count symbolic alts as non-ref
                            if(*alt == '<') any_alts = true;
                            // ignore empty / missing / symbolic alleles
                            continue;
                        }
                        any_alts = true;
                        ++allele_count;
                    }
                    else if(gt[g] == 0 && strict_homref)
                    {
                        ++allele_count;
                    }
                }
                if(vstart < previous_end && allele_count + previous_allele_count.get()[isample] > 2)
                {
                    if(all_warnings || !has_warned[WARNING::OVERLAP])
                    {
                        std::cerr << "[W] overlapping records at " << vchr << ":" << vstart << " for sample " << isample << "\n";
                    }
                    has_warned[WARNING::OVERLAP]++;
                }
                previous_allele_count.get()[isample] = allele_count;
            }

            if(any_alts)
            {
                ++n_nonref;
            }

            previous_end = vend;

            if(line->errcode)
            {
                error("Record at %s:%i will not translate into BCF. Check if the header is incomplete (error code %i)",
                      current_chr.c_str(), line->pos+1, line->errcode);
            }

            current_chr = vchr;

            if (message > 0 && (rcount % message) == 0)
            {
                std::cout << stringutil::formatPos(vchr.c_str(), line->pos) << "\n";
            }
            // count variants here
            ++rcount;
        }

        if(has_warned[WARNING::REFPADDING])
        {
            std::cerr << "[W] Variants that have >1 base of reference padding: " << has_warned[WARNING::REFPADDING] << "\n";
        }
        if(has_warned[WARNING::OVERLAP])
        {
            std::cerr << "[W] Variants that overlap on the reference allele: " << has_warned[WARNING::OVERLAP] << "\n";
        }

        std::cerr << "[I] Total VCF records:         " << rcount << "\n";
        std::cerr << "[I] Non-reference VCF records: " << n_nonref << "\n";
        bcf_sr_destroy(reader);
    }
    catch(std::runtime_error & e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch(std::logic_error & e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
