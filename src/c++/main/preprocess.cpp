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
 *  \brief VCF atomiser / preprocessor
 *
 *
 * \file xcmp.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"
#include "VariantInput.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"
#include "GraphReference.hh"
#include "DiploidCompare.hh"

#include "variant/VariantAlleleRemover.hh"
#include "variant/VariantAlleleSplitter.hh"
#include "variant/VariantTee.hh"
#include "variant/VariantAlleleNormalizer.hh"
#include "variant/VariantLocationAggregator.hh"
#include "variant/VariantHomrefSplitter.hh"
#include "variant/VariantAlleleUniq.hh"
#include "variant/VariantCallsOnly.hh"

#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>
#include <memory>

// error needs to come after boost headers.
#include "Error.hh"

using namespace variant;
using namespace haplotypes;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string ref_fasta;

    std::string chr = "";
    int64_t start = -1;
    int64_t end = -1;

    std::string file1;
    std::string sample1;

    std::string regions_bed = "";
    std::string targets_bed = "";

    std::string out_vcf = "";
    std::string out_errors = "";

    int64_t blimit = -1;
    bool progress = false;
    int progress_seconds = 10;

    bool preprocess = false;
    bool leftshift = false;
    bool haploid_X = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-vcf", po::value<std::string>(), "VCF files to preprocess (use file:sample for a specific sample column).")
            ("output-vcf,o", po::value<std::string>(), "Output variant comparison results to VCF.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("location,l", po::value<std::string>(), "The location to start at.")
            ("regions,R", po::value<std::string>(), "Use a bed file for getting a subset of regions (traversal via tabix).")
            ("targets,T", po::value<std::string>(), "Use a bed file for getting a subset of targets (streaming the whole file, ignoring things outside the bed regions).")
            ("progress", po::value<bool>(), "Set to true to output progress information.")
            ("haploid-x", po::value<bool>(), "Expand GTs on chrX: turn 1 into 1/1")
            ("progress-seconds", po::value<int>(), "Output progress information every n seconds.")
            ("limit", po::value<int64_t>(), "Maximum number of records to process.")
            ("preprocess-variants,V", po::value<bool>(), "Apply variant normalisations, trimming, realignment for complex variants (off by default).")
            ("leftshift,L", po::value<bool>(), "Left-shift indel alleles (off by default).")
        ;

        po::positional_options_description popts;
        popts.add("input-vcf", 1);

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
            std::cout << "preprocess version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-vcf"))
        {
            std::string vr = vm["input-vcf"].as< std::string >();
            std::vector<std::string> v;
            stringutil::split(vr, v, ":");
            // in case someone passes a ":"
            assert(v.size() > 0);

            file1 = v[0];
            sample1 = "";
            if(v.size() > 1)
            {
                sample1 = v[1];
            }
        }

        if (vm.count("output-vcf"))
        {
            out_vcf = vm["output-vcf"].as< std::string >();
        }

        if (vm.count("output-errors"))
        {
            out_errors = vm["output-errors"].as< std::string >();
        }

        if (vm.count("preprocess-variants"))
        {
            preprocess = vm["preprocess-variants"].as< bool >();
        }

        if (vm.count("leftshift"))
        {
            leftshift = vm["leftshift"].as< bool >();
        }

        if (vm.count("haploid-x"))
        {
            haploid_X = vm["haploid-x"].as< bool >();
        }

        if (vm.count("location"))
        {
            stringutil::parsePos(vm["location"].as< std::string >(), chr, start, end);
        }

        if (vm.count("regions"))
        {
            regions_bed = vm["regions"].as< std::string >();
        }

        if (vm.count("targets"))
        {
            targets_bed = vm["targets"].as< std::string >();
        }

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else if(preprocess)
        {
            error("Please specify a reference file name.");
        }

        if (vm.count("limit"))
        {
            blimit = vm["limit"].as< int64_t >();
        }

        if (vm.count("progress"))
        {
            progress = vm["progress"].as< bool >();
        }

        if (vm.count("progress-seconds"))
        {
            progress_seconds = vm["progress-seconds"].as< int >();
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
        VariantReader vr;
        vr.setReturnHomref(false);
        if(haploid_X)
        {
            vr.setFixChrXGTs(haploid_X);
        }

        if(regions_bed != "")
        {
            vr.setRegions(regions_bed.c_str(), true);
        }
        if(targets_bed != "")
        {
            vr.setTargets(targets_bed.c_str(), true);
        }

        int r1 = vr.addSample(file1.c_str(), sample1.c_str());

        vr.setApplyFilters(false, r1);

        VariantInput vi(
            ref_fasta.c_str(),
            preprocess || leftshift,          // bool leftshift
            true,          // bool refpadding
            true,                // bool trimalleles = false, (remove unused alleles)
            preprocess || leftshift,      // bool splitalleles = false,
            ( preprocess || leftshift ) ? 2 : 0,  // int mergebylocation = false,
            true,                // bool uniqalleles = false,
            true,                // bool calls_only = true,
            false,               // bool homref_split = false // this is handled by calls_only
            preprocess,          // bool primitives = false
            false,               // bool homref_output
            leftshift ? 1024 : 0, // int64_t leftshift_limit
            false
            );

        VariantProcessor & vp = vi.getProcessor();

        vp.setReader(vr, VariantBufferMode::buffer_block, 10*30);

        bool stop_after_chr_change = false;
        if(chr != "")
        {
            vp.rewind(chr.c_str(), start);
            stop_after_chr_change = true;
        }

        VariantWriter vw(out_vcf.c_str(), ref_fasta.c_str());
        vw.addHeader(vr);
        vw.setWriteFormats(true);
        std::list< std::pair<std::string, std::string> > files;
        vr.getSampleList(files);
        int sindex = 0;
        for(auto const & f : files)
        {
            if(f.second.empty())
            {
                vw.addSample((std::string("SAMPLE_") + std::to_string(sindex)).c_str());
            }
            else
            {
                vw.addSample(f.second.c_str());
            }
            ++sindex;
        }

        int64_t nrecs = 0;
        int64_t last_pos = std::numeric_limits<int64_t>::max();

        auto start_time = std::chrono::high_resolution_clock::now();
        auto last_time = std::chrono::high_resolution_clock::now();
        while(vp.advance())
        {
            if(blimit > 0 && nrecs++ > blimit)
            {
                // reached record limit
                break;
            }
            Variants & v = vp.current();

            if(end != -1 && (v.pos > end || (chr.size() != 0 && chr != v.chr)))
            {
                // reached end
                break;
            }

            if(stop_after_chr_change && chr.size() != 0 && chr != v.chr)
            {
                // reached end of chr
                break;
            }

            if(chr.size() == 0)
            {
                chr = v.chr;
            }

            chr = v.chr;
            vw.put(v);

            if(progress)
            {
                using namespace std;
                auto end_time = chrono::high_resolution_clock::now();
                auto secs = chrono::duration_cast<chrono::seconds>(end_time - last_time).count();

                if(secs > progress_seconds)
                {
                    auto secs_since_start = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
                    std::string mbps = "";
                    if(last_pos < v.pos)
                    {
                        mbps = " mpbs: ";
                        mbps += std::to_string(double(v.pos - last_pos) / double(secs_since_start) * 1e-6);
                    }
                    else
                    {
                        last_pos = v.pos;
                    }
                    last_time = end_time;

                    std::cerr << "[PROGRESS] Total time: " << secs_since_start << "s Pos: " << v.pos << mbps << "\n";
                }
            }
        }
    }
    catch(std::runtime_error &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch(std::logic_error &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

