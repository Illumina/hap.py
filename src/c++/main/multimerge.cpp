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
 *  \brief Merge many VCF/BCF files into one.
 *
 *  \details Mainly a test for VariantReader and VariantWriter
 *
 *
 * \file multimerge.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"

#include "VariantInput.hh"

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

    bool apply_filters = false;
    bool leftshift = false;
    bool trimalleles = false;
    bool splitalleles = false;
    int mergebylocation = false;
    bool uniqalleles = false;
    bool calls_only = true;
    bool homref_split = false;
    bool primitives = false;
    std::string homref_vcf = "";

    bool process_formats = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::vector< std::string> >(), "The input files")
            ("output-file,o", po::value<std::string>(), "The output file name.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("location,l", po::value<std::string>(), "Start location.")
            ("regions,R", po::value<std::string>(), "Use a bed file for getting a subset of regions (traversal via tabix).")
            ("targets,T", po::value<std::string>(), "Use a bed file for getting a subset of targets (streaming the whole file, ignoring things outside the bed regions).")
            ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
            ("message-every", po::value<int64_t>(), "Print a message every N records.")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("leftshift", po::value<bool>(), "Leftshift variant alleles.")
            ("trimalleles", po::value<bool>(), "Remove unused variant alleles.")
            ("splitalleles", po::value<bool>(), "Split and sort variant alleles.")
            ("merge-by-location", po::value<int>(), "Merge calls at the same location.")
            ("unique-alleles", po::value<bool>(), "Make alleles unique across a single line.")
            ("homref-split", po::value<bool>(), "Split homref blocks into per-nucleotide blocks.")
            ("homref-vcf-out", po::value<std::string>(), "Output split homref blocks as BCF/VCF.")
            ("calls-only", po::value<bool>(), "Remove homref blocks.")
            ("primitives", po::value<bool>(), "Split complex alleles into primitives via realignment.")
            ("process-split", po::value<bool>(), "Enables splitalleles, trimalleles, unique-alleles, leftshift.")
            ("process-full", po::value<bool>(), "Enables splitalleles, trimalleles, unique-alleles, leftshift, mergebylocation.")
            ("process-formats", po::value<bool>(), "Process GQ/DP/AD format fields.")
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

        if (vm.count("leftshift"))
        {
            leftshift = vm["leftshift"].as< bool >();
        }

        if (vm.count("trimalleles"))
        {
            trimalleles = vm["trimalleles"].as< bool >();
        }

        if (vm.count("splitalleles"))
        {
            splitalleles = vm["splitalleles"].as< bool >();
        }

        if (vm.count("merge-by-location"))
        {
            mergebylocation = vm["merge-by-location"].as< int >();
        }

        if (vm.count("unique-alleles"))
        {
            uniqalleles = vm["unique-alleles"].as< bool >();
        }

        if (vm.count("calls-only"))
        {
            calls_only = vm["calls-only"].as< bool >();
        }

        if (vm.count("homref-split"))
        {
            homref_split = vm["homref-split"].as< bool >();
        }

        if (vm.count("primitives"))
        {
            primitives = vm["primitives"].as< bool >();
        }

        if (vm.count("homref-vcf-out"))
        {
            homref_split = 1;
            homref_vcf = vm["homref-vcf-out"].as< std::string >();
        }

        if (vm.count("process-split"))
        {
            homref_split = true;
            trimalleles = true;
            splitalleles = true;
            uniqalleles = true;
            leftshift = true;
            calls_only = true;
        }

        if (vm.count("process-full"))
        {
            homref_split = true;
            trimalleles = true;
            splitalleles = true;
            uniqalleles = true;
            leftshift = true;
            calls_only = true;
            mergebylocation = 2;
            primitives = true;
        }

        if (vm.count("process-formats"))
        {
            process_formats = vm["process-formats"].as< bool >();
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
        std::shared_ptr<VariantWriter> p_homref_writer;
        if (homref_vcf.size() != 0)
        {
            p_homref_writer = std::make_shared<VariantWriter>(homref_vcf.c_str(), ref_fasta.c_str());
        }

        w.setWriteFormats(process_formats);
        if (p_homref_writer)
        {
            p_homref_writer->setWriteFormats(process_formats);
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

        std::list< std::pair<std::string, std::string> > samples;
        r.getSampleList(samples);

        std::set<std::string> samplenames;
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
            if (p_homref_writer)
            {
                p_homref_writer->addSample(sname.c_str());
            }
        }

        w.addHeader(r);
        if (p_homref_writer)
        {
            p_homref_writer->addHeader(r);
        }
        r.rewind(chr.c_str(), start);

        VariantInput vi(
            ref_fasta.c_str(),
            leftshift,           // bool leftshift
            true,                // bool refpadding
            trimalleles,         // bool trimalleles = false,
            splitalleles,        // bool splitalleles = false,
            mergebylocation,     // int mergebylocation = false,
            uniqalleles,         // bool uniqalleles = false,
            calls_only,          // bool calls_only = true,
            homref_split,        // bool homref_split = false
            primitives,          // bool primitives = false
            (bool)p_homref_writer// homref output
            );

        vi.getProcessor().setReader(r, VariantBufferMode::buffer_block, 100);

        VariantProcessor & proc = vi.getProcessor();
        VariantProcessor & proc_homref = vi.getProcessor(VariantInput::homref);

        int64_t rcount = 0;

        bool advance1 = proc.advance();
        bool advance2 = p_homref_writer ? proc_homref.advance() : false;

        while(advance1 || advance2)
        {
            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }

            if (advance1)
            {
                Variants & v = proc.current();
                if (chr.size() != 0 && chr != v.chr)
                {
                    // chromosome changed and location was given => abort
                    while(advance2)
                    {
                        Variants & v = proc_homref.current();
                        p_homref_writer->put(v);
                        advance2 = homref_split ? proc_homref.advance() : false;
                    }
                    break;
                }
                if(end != -1 && v.pos > end)
                {
                    break;
                }
                w.put(v);
                if(message > 0 && (rcount % message) == 0)
                {
                    std::cout << stringutil::formatPos(v.chr.c_str(), v.pos) << ": " << v << "\n";
                }
            }

            // make sure our homref variant output doesn't fill up all memory
            while(advance2 && p_homref_writer)
            {
                Variants & v = proc_homref.current();
                p_homref_writer->put(v);
                advance2 = proc_homref.advance();
            }

            advance1 = proc.advance();
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
