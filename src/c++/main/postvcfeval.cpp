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
 * Merge vcfeval output into a file that quantify understands
 *
 * \file postvcfeval.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Version.hh"
#include "Variant.hh"

#include "helpers/StringUtil.hh"

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <memory>

#include <htslib/hts.h>

#include "Error.hh"

#include "IntervalTree.h"

using namespace variant;

namespace std
{
    template<typename T, typename ...Args>
    unique_ptr<T> make_unique( Args&& ...args )
    {
        return unique_ptr<T>( new T( forward<Args>(args)... ) );
    }
}

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::string input_dir;
    std::string output_vcf;
    std::string ref_fasta;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-dir", po::value<std::string>(), "Path to a vcfeval output directory.")
            ("output-vcf", po::value<std::string>(), "Annotated VCF output file.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
        ;

        po::positional_options_description popts;
        popts.add("input-dir", 1);
        popts.add("output-vcf", 1);

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
            std::cout << "postvcfeval version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-dir"))
        {
            input_dir = vm["input-dir"].as<std::string>();
        }
        else
        {
            error("Please specify an input directory");
        }

        if (vm.count("output-vcf"))
        {
            output_vcf = vm["output-vcf"].as< std::string >();
        }

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else
        {
            error("To write an output VCF, you need to specify a reference file, too.");
        }

        if (output_vcf == "")
        {
            std::cerr << "Please specify an output file.\n";
            return 1;
        }
    }
    catch (po::error & e)
    {
        std::cerr << e.what() << "\n";
        return 1;
    }

    try
    {
        VariantReader r;
        r.setApplyFilters(false);

        boost::filesystem::path p(input_dir);
        int in_ix_fp  = r.addSample((p / "fp.vcf.gz").c_str(), "");
        int in_ix_fn  = r.addSample((p / "fn.vcf.gz").c_str(), "");
        int in_ix_tp  = r.addSample((p / "tp.vcf.gz").c_str(), "");
        int in_ix_tpb = r.addSample((p / "tp-baseline.vcf.gz").c_str(), "");

        std::unique_ptr<VariantWriter> writer = std::make_unique<VariantWriter>(output_vcf.c_str(), ref_fasta.c_str());
        writer->addSample("TRUTH");
        writer->addSample("QUERY");
        writer->addHeader(r);
        writer->addHeader("##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
        writer->addHeader("##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");

        int64_t rcount = 0;

        while(r.advance())
        {
            Variants & v = r.current();

            if(!v.calls[in_ix_fp].isNocall())
            {
                Variants out_vars;
                out_vars = v;
                out_vars.calls.clear();
                out_vars.calls.resize(2);

                out_vars.calls[1] = v.calls[in_ix_fp];
                if(!out_vars.info.empty()) out_vars.info += ";";
                out_vars.info += "type=FP";
                out_vars.info += ";kind=missing";
                writer->put(out_vars);
            }

            if(!v.calls[in_ix_fn].isNocall())
            {
                Variants out_vars;
                out_vars = v;
                out_vars.calls.clear();
                out_vars.calls.resize(2);

                out_vars.calls[0] = v.calls[in_ix_fn];
                if(!out_vars.info.empty()) out_vars.info += ";";
                out_vars.info += "type=FN";
                out_vars.info += ";kind=missing";
                writer->put(out_vars);
            }

            if(!v.calls[in_ix_tp].isNocall() || !v.calls[in_ix_tpb].isNocall())
            {
                Variants out_vars;
                out_vars = v;
                out_vars.calls.clear();
                out_vars.calls.resize(2);

                out_vars.calls[0] = v.calls[in_ix_tp];
                out_vars.calls[1] = v.calls[in_ix_tpb];
                if(!out_vars.info.empty()) out_vars.info += ";";
                out_vars.info += "type=TP";
                out_vars.info += ";kind=vcfeval";
                writer->put(out_vars);
            }
            ++rcount;
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

