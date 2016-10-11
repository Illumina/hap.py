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
 * \brief Split a set of VCFs into blocks with variants not closer than a given window length
 *
 * \file blocksplit.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"

#include <fstream>
#include <limits>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

// error needs to come after program_options.
#include "Error.hh"

using namespace variant;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::vector<std::string> files, samples;

    std::string regions_bed = "";
    std::string targets_bed = "";
    std::string out_bed = "";

    // limits
    std::string chr = "";
    int64_t start = -1;
    int64_t end = -1;

    int64_t rlimit = -1;
    int64_t message = -1;
    int64_t window = 30;
    bool apply_filters = false;

    int nblocks = 32;
    int nvars = 100;

    bool verbose = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::vector<std::string> >(), "The input VCF/BCF file(s) (use file:sample to specify a sample)")
            ("output,o", po::value<std::string>(), "Write a bed file giving the locations of overlapping blocks (use - for stdout).")
            ("regions,R", po::value<std::string>(), "Use a bed file for getting a subset of regions (traversal via tabix).")
            ("targets,T", po::value<std::string>(), "Use a bed file for getting a subset of targets (streaming the whole file, ignoring things outside the bed regions).")
            ("location,l", po::value<std::string>(), "The location / subset.")
            ("limit-records,L", po::value<int64_t>(), "Maximum number of records to process")
            ("message-every,m", po::value<int64_t>(), "Print a message every N records.")
            ("window,w", po::value<int64_t>(), "Overlap window length.")
            ("nblocks,b", po::value<int>(), "Maximum number of blocks to break into (32).")
            ("nvars,v", po::value<int>(), "Minimum number of variants per block (100).")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("verbose", po::value<bool>(), "Verbose output.")
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
            std::cout << "blocksplit version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-file"))
        {
            std::vector<std::string> fs = vm["input-file"].as< std::vector<std::string> >();

            for(std::string const & s : fs)
            {
                std::vector<std::string> v;
                stringutil::split(s, v, ":");
                std::string filename, sample = "";

                // in case someone passes a ":"
                assert(v.size() > 0);

                filename = v[0];

                if(v.size() > 1)
                {
                    sample = v[1];
                }

                files.push_back(filename);
                samples.push_back(sample);
            }
        }

        if(files.size() == 0)
        {
            error("Please specify at least one input file.");
        }

        if (vm.count("output"))
        {
            out_bed = vm["output"].as< std::string >();
        }
        else
        {
            out_bed = "-";
        }

        if (vm.count("regions"))
        {
            regions_bed = vm["regions"].as< std::string >();
        }

        if (vm.count("targets"))
        {
            targets_bed = vm["targets"].as< std::string >();
        }

        if (vm.count("verbose"))
        {
            verbose = vm["verbose"].as<bool>();
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

        if (vm.count("window"))
        {
            window = vm["window"].as< int64_t >();
        }

        if (vm.count("nblocks"))
        {
            nblocks = vm["nblocks"].as< int >();
        }

        if (vm.count("nvars"))
        {
            nvars = vm["nvars"].as< int >();
        }

    }
    catch (po::error & e)
    {
        std::cerr << e.what() << "\n";
        return 1;
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

    try
    {
        bcf_srs_t * reader = bcf_sr_init();
        reader->require_index = 1;
        reader->collapse = COLLAPSE_NONE;
        reader->streaming = 0;
        if(!regions_bed.empty())
        {
            int result = bcf_sr_set_regions(reader, regions_bed.c_str(), 1);
            if(result < 0)
            {
                error("Failed to set regions string %s.", regions_bed.c_str());
            }
        }
        if(!targets_bed.empty())
        {
            int result = bcf_sr_set_targets(reader, targets_bed.c_str(), 1, 1);
            if(result < 0)
            {
                error("Failed to set targets string %s.", targets_bed.c_str());
            }
        }
        for(auto const & file : files)
        {
            if (!bcf_sr_add_reader(reader, file.c_str()))
            {
                error("Failed to open or file not indexed: %s\n", file.c_str());
            }
        }

        bool stop_after_chr_change = false;
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
            if(success == -reader->nreaders)
            {
                // cannot seek -> return no output
                // write blocks
                if(out_bed == "-" || out_bed == "")
                {
                    return 0;
                }
                else
                {
                    if(verbose)
                    {
                        std::cerr << "Writing to " << out_bed << "\n";
                    }
                    volatile std::ofstream f(out_bed.c_str());
                    return 0;
                }
            }
            stop_after_chr_change = true;
        }

        int64_t rcount = 0;
        int64_t last_end = -1;
        int64_t vars = 0, total_vars = 0;

        struct Breakpoint
        {
            std::string chr;
            int64_t pos;
            int64_t vars;
        };

        std::list< Breakpoint > breakpoints;

        const auto add_bp = [&breakpoints, nvars, nblocks, &chr, &vars, verbose](int64_t bp)
        {
            if (vars > nvars)
            {
                if(verbose)
                {
                    std::cerr << "Break point at " << chr << ":" << bp << " (" << vars << " variants)" << "\n";
                }
                breakpoints.push_back(Breakpoint{chr, bp, vars});
                vars = 0;
            }
        };

        std::string firstchr;

        int nl = 1;
        while(nl)
        {
            nl = bcf_sr_next_line(reader);
            if (nl <= 0)
            {
                break;
            }

            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }

            std::string v_chr;
            int64_t v_pos = -1, v_end = -1;
            for(int isample = 0; isample < reader->nreaders; ++isample)
            {
                if(!bcf_sr_has_line(reader, isample))
                {
                    continue;
                }
                bcf_hdr_t * hdr = reader->readers[isample].header;
                bcf1_t * line = reader->readers[isample].buffer[0];
                bcf_unpack(line, BCF_UN_FLT);

                if(apply_filters)
                {
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

                bool call_this_pos = false;
                for(int rsample = 0; rsample < line->n_sample; ++rsample)
                {
                    int gts[MAX_GT];
                    int ngt = 0;
                    bool phased = false;
                    bcfhelpers::getGT(hdr, line, rsample, gts, ngt, phased);
                    for(int j = 0; j < ngt; ++j)
                    {
                        if(gts[j] > 0)
                        {
                            call_this_pos = true;
                            break;
                        }
                    }
                    if(call_this_pos)
                    {
                        break;
                    }
                }

                if(!call_this_pos)
                {
                    continue;
                }

                v_chr = bcfhelpers::getChrom(hdr, line);
                // rely on synced_reader to give us records on the same chr
                int64_t this_v_pos = -1;
                int64_t this_v_end = -1;
                try
                {
                    bcfhelpers::getLocation(hdr, line, this_v_pos, this_v_end);
                }
                catch(bcfhelpers::importexception const & e)
                {
                    std::cerr << e.what() << "\n";
                    continue;
                }
                if(v_pos < 0)
                {
                    v_pos = this_v_pos;
                }
                else
                {
                    v_pos = std::min(this_v_pos, v_pos);
                }
                v_end = std::max(this_v_end, v_end);
            }

            if(v_chr.empty() || v_pos < 0 || v_end < 0)
            {
                continue;
            }

            if(end != -1 && ( (v_pos > end) || (chr != "" && v_chr != chr)) )
            {
                break;
            }
            if(stop_after_chr_change && chr != "" && v_chr != chr)
            {
                break;
            }
            if (firstchr.size() == 0)
            {
                firstchr = v_chr;
            }
            if(chr != "" && v_chr != chr)
            {
                last_end = -1;
            }
            chr = v_chr;

            if(message > 0 && (rcount % message) == 0)
            {
                std::cerr << "From " << chr << ":" << last_end << " ("
                          << breakpoints.size() << " bps, " << vars << " vars)"
                          << " -- " << v_chr << ":" << v_pos << "-" << v_end << "\n";
            }

            vars++;
            total_vars++;

            if(last_end >= 0 && v_pos > last_end + window) // can split here
            {
                add_bp(last_end);
            }
            last_end = std::max(last_end, v_end);

            ++rcount;
        }

        // write blocks
        std::ostream * outputfile = NULL;

        if(out_bed == "-" || out_bed == "")
        {
            outputfile = &std::cout;
        }
        else
        {
            if(verbose)
            {
                std::cerr << "Writing to " << out_bed << "\n";
            }
            outputfile = new std::ofstream(out_bed.c_str());
        }

        chr = firstchr;
        start = 0;
        int64_t vpb = 0;
        int64_t target_vpb = std::max(nvars, ((int)total_vars) / nblocks);

        if(end <= 0)
        {
            end = std::numeric_limits<int>::max();
        }

        for (auto & b : breakpoints)
        {
            if (chr != b.chr)
            {
                *outputfile << chr << "\t" << start << "\t" << std::max(start + window + 1, end) << "\n";
                chr = b.chr;
                start = 1;
                vpb = 0;
            }
            vpb += b.vars;
            if(vpb > target_vpb)
            {
                *outputfile << chr << "\t" << start << "\t" << b.pos + window + 1 << "\n";
                start = b.pos + window + 1;
                vpb = 0;
            }
        }
        if(chr != "")
        {
            *outputfile << chr << "\t" << start << "\t" << std::max(start + window + 1, end) << "\n";
        }

        if(out_bed != "-" && out_bed != "")
        {
            delete outputfile;
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

