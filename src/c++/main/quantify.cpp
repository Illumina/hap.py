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
 * Count variants in a VCF file
 *
 * \file quantify.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
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

#include "Error.hh"

#include "BlockQuantify.hh"

using namespace variant;


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::vector<std::string> files;
    std::string output;
    std::string output_vcf;
    std::string ref_fasta;
    std::string only_regions;

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;
    int64_t rlimit = -1;

    int64_t message = -1;

    bool apply_filters = false;
    bool count_homref = false;
    bool output_vtc = false;

    QuantifyRegions regions;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::vector< std::string> >(), "The input files")
            ("output-file,o", po::value<std::string>(), "The output file name (JSON Format).")
            ("output-vcf,v", po::value<std::string>(), "Annotated VCF file (with bed annotations).")
            ("reference,r", po::value<std::string>(), "The reference fasta file (needed only for VCF output).")
            ("location,l", po::value<std::string>(), "Start location.")
            ("regions,R", po::value< std::vector<std::string> >(),
                "Region bed file. You can attach a label by prefixing with a colon, e.g. -R FP2:false-positives-type2.bed")
            ("only,O", po::value< std::string >(), "Bed file of locations (equivalent to -R in bcftools)")
            ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
            ("message-every", po::value<int64_t>(), "Print a message every N records.")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("count-homref", po::value<bool>(), "Count homref locations.")
            ("output-vtc", po::value<bool>(), "Output variant types counted (debugging).")
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
            std::cout << "quantify version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-file"))
        {
            files = vm["input-file"].as< std::vector<std::string> >();
        }

        if (vm.count("output-file"))
        {
            output = vm["output-file"].as< std::string >();
        }

        if (vm.count("output-vcf"))
        {
            output_vcf = vm["output-vcf"].as< std::string >();
        }

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else if(output_vcf != "")
        {
            error("To write an output VCF, you need to specify a reference file, too.");
        }

        if (vm.count("location"))
        {
            stringutil::parsePos(vm["location"].as< std::string >(), chr, start, end);
        }

        if (vm.count("only"))
        {
            only_regions = vm["only"].as< std::string >();
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

        if (vm.count("count-homref"))
        {
            count_homref = vm["count-homref"].as< bool >();
        }

        if (vm.count("output-vtc"))
        {
            output_vtc = vm["output-vtc"].as< bool >();
        }

        if(files.size() == 0)
        {
            std::cerr << "Please specify at least one input file / sample.\n";
            return 1;
        }

        if (output == "")
        {
            std::cerr << "Please specify an output file.\n";
            return 1;
        }

        if (vm.count("regions"))
        {
            std::vector<std::string> rnames = vm["regions"].as< std::vector<std::string> >();
            regions.load(rnames);
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
        if(!only_regions.empty()) {
            r.setRegions(only_regions.c_str(), true);
        }
        r.setApplyFilters(apply_filters);

        for(std::string const & f : files)
        {
            std::vector<std::string> v;
            stringutil::split(f, v, ":");

            std::string filename, sample = "";

            // in case someone passes a ":"
            assert(v.size() > 0);

            filename = v[0];

            if(v.size() > 1)
            {
                sample = v[1];
            }
            std::cerr << "Adding file '" << filename << "' / sample '" << sample << "'" << "\n";
            r.addSample(filename.c_str(), sample.c_str());
        }

        r.rewind(chr.c_str(), start);

        std::shared_ptr<VariantWriter> writer;

        if (output_vcf != "")
        {
            writer = std::make_shared<VariantWriter>(output_vcf.c_str(), ref_fasta.c_str());
            std::set<std::string> samplenames;
            std::list< std::pair<std::string, std::string> > samples;
            r.getSampleList(samples);
            for (auto const & p : samples)
            {
                std::string sname = p.second;
                if (sname == "")
                {
                    sname = boost::filesystem::path(p.first).stem().string();
                }
                int i = 1;
                while (samplenames.count(sname))
                {
                    sname = p.second + "." + std::to_string(i++);
                }
                samplenames.insert(sname);
                std::cerr << "Writing '" << p.first << ":" << p.second << "' as sample '" << sname << "'" << "\n";
                writer->addSample(sname.c_str());
            }
            writer->addHeader(r);
            writer->addHeader("##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">");
            writer->addHeader("##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">");
            writer->addHeader("##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
            writer->addHeader("##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
            writer->addHeader("##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">");
            if(output_vtc) {
                writer->addHeader("##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">");
            }
            writer->addHeader("##INFO=<ID=T_VT,Number=1,Type=String,Description=\"High-level variant type in truth (SNP|INDEL).\">");
            writer->addHeader("##INFO=<ID=Q_VT,Number=1,Type=String,Description=\"High-level variant type in query (SNP|INDEL).\">");
            writer->addHeader("##INFO=<ID=T_LT,Number=1,Type=String,Description=\"High-level location type in truth (het|hom|hetalt).\">");
            writer->addHeader("##INFO=<ID=Q_LT,Number=1,Type=String,Description=\"High-level location type in query (het|hom|hetalt).\">");
        }

        /** local function to count variants in all samples */
        int64_t rcount = 0;
        std::string current_chr = "";
        BlockQuantify bq(r, regions, ref_fasta, output_vtc, count_homref);
        while(r.advance())
        {
            Variants & v = r.current();

            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }

            if(end != -1 && ((!current_chr.empty() && chr != current_chr) || v.pos > end))
            {
                break;
            }
            current_chr = v.chr;
            bq.add(v);
            if (message > 0 && (rcount % message) == 0) {
                std::cout << stringutil::formatPos(v.chr.c_str(), v.pos) << ": " << v << "\n";
            }
            // count variants here
            ++rcount;
        }

        bq.count();

        if(writer)
        {
            auto const & variants = bq.getVariants();
            for(auto & v : variants)
            {
                writer->put(v);
            }
        }

        auto const & count_map = bq.getCounts();

        Json::Value counts;
        for (auto & p : count_map)
        {
            counts[p.first.c_str()] = p.second.write();
        }

        Json::FastWriter fastWriter;
        if (output == "" || output == "-")
        {
            std::cout << fastWriter.write(counts);
        }
        else
        {
            std::ofstream o(output);
            o << fastWriter.write(counts);
        }

        // also print / write
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
