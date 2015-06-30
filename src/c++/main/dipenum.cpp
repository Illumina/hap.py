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
 * \brief Enumerate possible diploid Haplotypes from a VCF
 *
 * \file dipenum.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"
#include "Haplotype.hh"
#include "DiploidReference.hh"
#include "VariantInput.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"

#include <fstream>

// error needs to come after program_options.
#include "Error.hh"

using namespace variant;
using namespace haplotypes;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string ref_fasta;
    std::string file;
    std::string sample;

    std::string out_fasta = "";

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;

    // = max 12 unphased hets in segment
    int max_n_haplotypes = 4096;

    bool apply_filters = true;
    bool preprocess = true;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::string>(), "The input VCF/BCF file (use file:sample to specify a sample)")
            ("output", po::value<std::string>(), "Write a file with all possible haplotypes.")
            ("location,l", po::value<std::string>(), "The location / subset.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("max-n-haplotypes", po::value<int>(), "Maximum number of haplotypes to enumerate.")
            ("apply-filters,f", po::value<int>(), "Apply filters in VCF (default to 1)")
            ("preprocess,P", po::value<bool>(), "Preprocess variants")
        ;

        po::positional_options_description popts;
        popts.add("input-file", 1);

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
            std::cout << "dipenum version " << HAPLOTYPES_VERSION << "\n";
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
            std::vector<std::string> v;
            stringutil::split(file, v, ":");

            // in case someone passes a ":"
            assert(v.size() > 0);

            file = v[0];

            if(v.size() > 1)
            {
                sample = v[1];
            }
        }

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else
        {
            error("Please specify a reference file name.");
        }

        if (vm.count("output"))
        {
            out_fasta = vm["output"].as< std::string >();
        }

        if (vm.count("location"))
        {
            stringutil::parsePos(vm["location"].as< std::string >(), chr, start, end);
        }
        else
        {
            error("Please specify a location.");
        }

        if (vm.count("max-n-haplotypes"))
        {
            max_n_haplotypes = vm["max-n-haplotypes"].as< int >();
        }

        if (vm.count("apply-filters"))
        {
            apply_filters = vm["apply-filters"].as< int >() != 0;
        }

        if (vm.count("preprocess"))
        {
            preprocess = vm["preprocess"].as<bool>() != 0;
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
        VariantReader r;
        r.setApplyFilters(apply_filters);
        int ix = r.addSample(file.c_str(), sample.c_str());
        VariantInput vi(
            ref_fasta.c_str(),
            preprocess,           // bool leftshift
            false,           // bool refpadding
            false,          // bool trimalleles = false,
            false,          // bool splitalleles = false,
            0,              // int mergebylocation = false,
            preprocess,           // bool uniqalleles = false,
            false,          // bool calls_only = true,
            false           // bool homref_split = false
        );
        vi.getProcessor().setReader(r, VariantBufferMode::buffer_block, 100);

        DiploidReference dr(ref_fasta.c_str());
        dr.setNPaths(max_n_haplotypes);

        std::ostream * out;
        if(out_fasta != "" && out_fasta != "-")
        {
            out = new std::ofstream(out_fasta.c_str());
        }
        else
        {
            out = &std::cout;
        }

        std::list<Variants> vars;
        vi.get(chr.c_str(), start, end, vars);
        dr.setRegion(chr.c_str(), start, end, vars, ix);

        while(dr.hasNext())
        {
            DiploidRef & hp (dr.next());

            *out << hp << "\n";

            dr.advance();
        }


        if(out_fasta != "" && out_fasta != "-")
        {
            delete out;
        }
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

