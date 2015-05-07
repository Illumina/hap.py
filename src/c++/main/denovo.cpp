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
 * \brief Find denovo variants in one sample vs. many.
 *
 * \file denovo.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */


#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"

#include "variant/VariantAlleleRemover.hh"
#include "variant/VariantAlleleSplitter.hh"
#include "variant/VariantAlleleNormalizer.hh"
#include "variant/VariantLocationAggregator.hh"
#include "variant/VariantHomrefSplitter.hh"
#include "variant/VariantAlleleUniq.hh"
#include "variant/VariantCallsOnly.hh"

#include "denovo/MarkDenovo.hh"

#include "helpers/StringUtil.hh"

#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <fstream>
#include <set>

#include "Error.hh"

using namespace variant;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::vector<std::string> files;
    std::string output;
    std::string ref_fasta;
    std::string regions_bed = "";
    std::string targets_bed = "";

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;
    int64_t rlimit = -1;

    int64_t message = -1;

    bool apply_filters = true;
    bool output_all = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::vector< std::string> >(), "The input files / samples. First sample is the child / sample to find denovo variants in. "
             "The next two are the parents. All the following ones are used to remove variants.")
            ("output-file,o", po::value<std::string>(), "The output file name.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("location,l", po::value<std::string>(), "Start location.")
            ("regions,R", po::value<std::string>(), "Use a bed file for getting a subset of regions (traversal via tabix).")
            ("targets,T", po::value<std::string>(), "Use a bed file for getting a subset of targets (streaming the whole file, ignoring things outside the bed regions).")
            ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
            ("message-every", po::value<int64_t>(), "Print a message every N records.")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("output-all,A", po::value<bool>(), "Output all variants.")
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
            std::cout << "multimerge version " << HAPLOTYPES_VERSION << "\n";
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

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else
        {
            error("Please specify a reference file name.");
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

        if (vm.count("output-all"))
        {
            output_all = vm["output-all"].as< bool >();
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
    }
    catch (po::error & e)
    {
        std::cerr << e.what() << "\n";
        return 1;
    }

    try
    {
        VariantReader r;
        if(regions_bed != "")
        {
            r.setRegions(regions_bed.c_str(), true);
        }
        if(targets_bed != "")
        {
            r.setTargets(targets_bed.c_str(), true);
        }

        VariantWriter w(output.c_str(), ref_fasta.c_str());
        w.setWriteFormats(true);
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

        std::list< std::pair<std::string, std::string> > samples;
        r.getSampleList(samples);

        if(samples.size() < 3)
        {
            error("I need at least three samples.");
        }

        std::set<std::string> samplenames;
        int counter = 0;
        for (auto const & p : samples)
        {
            std::string sname = p.second;
            if (sname == "")
            {
                sname = boost::filesystem::path(p.first).stem().string();
            }
            switch(counter)
            {
                case 0:
                    sname = std::string("CHILD.") + sname;
                    break;
                case 1:
                case 2:
                    sname = std::string("PARENT.") + std::to_string(counter) + "." + sname;
                    break;
                default:
                    sname = std::string("SIBLING.") + std::to_string(counter-2) + "." + sname;
                    break;
            }
            int i = 1;
            while (samplenames.count(sname))
            {
                if(p.second == "")
                {
                    sname = boost::filesystem::path(p.first).stem().string() + "." + std::to_string(i++);
                }
                else
                {
                    sname = p.second + "." + std::to_string(i++);
                }
            }
            samplenames.insert(sname);
            std::cerr << "Writing '" << p.first << ":" << p.second << "' as sample '" << sname << "'" << "\n";
            w.addSample(sname.c_str());
            ++counter;
        }

        w.addHeader(r);
        w.addHeader("##INFO=<ID=denovo,Number=.,Type=String,Description=\"Type of denovo variant = child_only|parents_only|child_partial\">");

        bool stop_after_chr_change = false;
        if(chr.size() != 0)
        {
            r.rewind(chr.c_str(), start);
            stop_after_chr_change = true;
        }

        VariantProcessor proc;

        // split homref blocks
        VariantHomrefSplitter homref_splitter;
        proc.addStep(homref_splitter);

        VariantCallsOnly calls_only;
        proc.addStep(calls_only);

        VariantAlleleRemover allele_remover;
        proc.addStep(allele_remover);

        VariantAlleleSplitter allele_splitter;
        proc.addStep(allele_splitter);

        VariantAlleleNormalizer allele_normalizer;
        allele_normalizer.setReference(ref_fasta.c_str());
        allele_normalizer.setEnableRefPadding(true);
        allele_normalizer.setEnableHomrefVariants(true);
        proc.addStep(allele_normalizer);

        VariantLocationAggregator merger;
        merger.setAggregationType(VariantLocationAggregator::aggregate_hetalt);
        proc.addStep(merger);

        VariantAlleleUniq allele_uniq;
        proc.addStep(allele_uniq);

        MarkDenovo mdn;
        proc.addStep(mdn);

        proc.setReader(r, VariantBufferMode::buffer_block, 100);

        int64_t rcount = 0;

        while(proc.advance())
        {
            Variants & v = proc.current();

            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }

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

            if(output_all || v.info.find_first_of("denovo") != std::string::npos)
            {
                w.put(v);
            }

            if(message > 0 && (rcount % message) == 0)
            {
                std::cout << stringutil::formatPos(v.chr.c_str(), v.pos) << ": " << v << "\n";
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
