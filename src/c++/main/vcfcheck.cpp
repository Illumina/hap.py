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

// error needs to come after boost headers.
#include "Error.hh"

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
                ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
                ("message-every", po::value<int64_t>(), "Print a message every N records.")
                ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
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
                ++vstart;
                if(refpadding > 1)
                {
                    std::cerr << "[W] variant at " << vchr << ":" << vstart << " has more than one base of reference padding \n";
                }
            }

            if(vstart < previous_end)
            {
                std::cerr << "[W] overlapping records at " << vchr << ":" << vstart << "\n";
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
