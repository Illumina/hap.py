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
 * \brief Haplotype comparison of VCF regions
 *
 * \file hapcmp.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"
#include "Haplotype.hh"
#include "GraphReference.hh"
#include "DiploidCompare.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"

#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>

// error needs to come after program_options. 
#include "Error.hh"

using namespace variant;
using namespace haplotypes;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string ref_fasta;

    std::string regions;

    std::string file1;
    std::string sample1;
    std::string file2;
    std::string sample2;

    std::string out_bed = "";
    std::string out_errors = "";
    std::string out_diffs = "";

    // = max 12 unphased hets in segment
    int max_n_haplotypes = 4096;
    int64_t blimit = -1;
    bool progress = false;
    int progress_seconds = 10;

    bool output_sequences = false;
    bool apply_filters = false;
    bool do_alignment = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-regions", po::value<std::string>(), "The input bed file specifying haplotype block regions (use - for stdin).")
            ("input-vcfs", po::value<std::vector<std::string> >(), "Two VCF files to compare (use file:sample for a specific sample column).")
            ("output-bed,b", po::value<std::string>(), "Output block results as bed files (default is to output to stdout).")
            ("output-errors,e", po::value<std::string>(), "Output failure information.")            
            ("output-diffs,d", po::value<std::string>(), "Output shared and different variants to a mJSON file (one json record per line, default is to not output diffs).")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("max-n-haplotypes,n", po::value<int>(), "Maximum number of haplotypes to enumerate.")
            ("output-sequences", po::value<bool>(), "Set to true to output haplotype sequences.")
            ("progress", po::value<bool>(), "Set to true to output progress information.")
            ("progress-seconds", po::value<int>(), "Output progress information every n seconds.")
            ("limit,l", po::value<int64_t>(), "Maximum number of haplotype blocks to process.")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("do-alignment", po::value<bool>(), "Perform alignments on mismatching haplotypes to find best approximate match.")
        ;

        po::positional_options_description popts;
        popts.add("input-regions", 1);
        popts.add("input-vcfs", 2);

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
            std::cout 
                << "hapcmp version " << HAPLOTYPES_VERSION 
                << "\n";
            return 0;
        }

        if (vm.count("help")) 
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-regions"))
        {
            regions = vm["input-regions"].as< std::string >();
        }
        else
        {
            error("Please specify input regions.");
        }

        if (vm.count("input-vcfs")) 
        {
            std::vector<std::string> vr = vm["input-vcfs"].as< std::vector<std::string> >();

            if(vr.size() != 2)
            {
                error("Please pass exactly two vcf file names for comparison.");
            }

            std::vector<std::string> v;
            stringutil::split(vr[0], v, ":");
            // in case someone passes a ":"
            assert(v.size() > 0);

            file1 = v[0];
            sample1 = "";
            if(v.size() > 1)
            {
                sample1 = v[1];
            }

            v.clear();
            stringutil::split(vr[1], v, ":");
            // in case someone passes a ":"
            assert(v.size() > 0);

            file2 = v[0];
            sample2 = "";
            if(v.size() > 1)
            {
                sample2 = v[1];
            }
        }

        if (vm.count("output-bed"))
        {
            out_bed = vm["output-bed"].as< std::string >();
        }

        if (vm.count("output-diffs"))
        {
            out_diffs = vm["output-diffs"].as< std::string >();
        }

        if (vm.count("output-errors"))
        {
            out_errors = vm["output-errors"].as< std::string >();
        }

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else
        {
            error("Please specify a reference file name.");
        }

        if (vm.count("max-n-haplotypes"))
        {
            max_n_haplotypes = vm["max-n-haplotypes"].as< int >();
        }

        if (vm.count("limit"))
        {
            blimit = vm["limit"].as< int64_t >();
        }

        if (vm.count("output-sequences"))
        {
            output_sequences = vm["output-sequences"].as< bool >();
        }

        if (vm.count("progress"))
        {
            progress = vm["progress"].as< bool >();
        }

        if (vm.count("progress-seconds"))
        {
            progress_seconds = vm["progress-seconds"].as< int >();
        }

        if (vm.count("apply-filters"))
        {
            apply_filters = vm["apply-filters"].as< bool >();
        }

        if (vm.count("do-alignment"))
        {
            do_alignment = vm["do-alignment"].as< bool >();
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
        vr.setApplyFilters(apply_filters);
        int ix1 = vr.addSample(file1.c_str(), sample1.c_str());
        int ix2 = vr.addSample(file2.c_str(), sample2.c_str());

        DiploidCompare dc(ref_fasta.c_str());
        dc.setMaxHapEnum(max_n_haplotypes);
        dc.setDoAlignments(do_alignment);

        std::istream * in = NULL;
        if(regions != "-")
        {
            in = new std::ifstream(regions.c_str());
        }
        else
        {
            in = &std::cin;
        }

        std::ostream * bed_out_stream = NULL;
        if(out_bed == "" || out_bed == "-")
        {
            bed_out_stream = &std::cout;
        }
        else
        {
            bed_out_stream = new std::ofstream(out_bed.c_str());
        }

        std::ostream * error_out_stream = NULL;
        if(out_errors == "" || out_errors == "-")
        {
            error_out_stream = &std::cerr;
        }
        else
        {
            error_out_stream = new std::ofstream(out_errors.c_str());            
        }

        std::ostream * diff_out_stream = NULL;
        if(out_diffs != "")
        {
            diff_out_stream = new std::ofstream(out_diffs.c_str());
        }

        Json::FastWriter fastWriter;

        int64_t nhb = 0;
        int64_t last_pos = std::numeric_limits<int64_t>::max();

        auto start_time = std::chrono::high_resolution_clock::now();
        auto last_time = std::chrono::high_resolution_clock::now();

        while(in->good())
        {
            if(blimit > 0 && nhb++ > blimit)
            {
                // reached record limit
                break;
            }

            std::string line;
            std::getline(*in, line);

            std::vector<std::string> v;

            stringutil::split(line, v, "\t");
            // we want >= 3 columns
            if(v.size() > 3)
            {
                std::string chr = v[0];
                int64_t start = std::stoll(v[1]);
                int64_t end = std::stoll(v[2]) - 1;

                if(progress)
                {
                    using namespace std;
                    auto end_time = chrono::high_resolution_clock::now();
                    auto secs = chrono::duration_cast<chrono::seconds>(end_time - last_time).count();

                    if(secs > progress_seconds)
                    {
                        auto secs_since_start = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
                        std::string mbps = "";
                        if(last_pos < end)
                        {
                            mbps = " mpbs: ";
                            mbps += std::to_string(double(end - last_pos) / double(secs_since_start) * 1e-6);
                        }
                        else
                        {
                            last_pos = end;                            
                        }
                        last_time = end_time;

                        std::cerr << "[PROGRESS] Total time: " << secs_since_start << "s Pos: " << end << mbps << "\n";
                    }                    
                }

                try
                {
                    std::ostringstream partial_bed;
                    vr.rewind(chr.c_str(), start);
                    std::list<Variants> vars;

                    while(vr.advance())
                    {
                        Variants & v = vr.current();
                        if(chr != v.chr)
                        {
                            // break on change of chr
                            break;
                        }
                        if(v.pos > end)
                        {
                            break;
                        }

                        vars.push_back(v);
                    }
                    dc.setRegion(chr.c_str(), start, end, vars, ix1, ix2);
                    DiploidComparisonResult const & dcr = dc.getResult();

                    printDiploidComparisonResult(partial_bed, dcr, output_sequences);
                    partial_bed << "\n";

                    if (diff_out_stream)
                    {
                        *diff_out_stream << fastWriter.write(toJson(dcr));
                    }

                    *bed_out_stream << partial_bed.str();
                }
                catch(std::runtime_error &e)
                {
                    *bed_out_stream << chr << "\t" << start << "\t" << end+1 << "\t" << dco_unknown << "\t.\t.\t.\t.\t.\t.\t.\n";
                    *error_out_stream << "[E] Error comparing block at " << chr << ":" << start << "-" << end << " - " << e.what() << "\n";
                }
                catch(std::logic_error &e)
                {
                    *bed_out_stream << chr << "\t" << start << "\t" << end+1 << "\t" << dco_unknown << "\t.\t.\t.\t.\t.\t.\t.\n";
                    *error_out_stream << "[E] Logic error comparing block at " << chr << ":" << start << "-" << end << " - " << e.what()  << "\n";
                }

            }
            else
            {
                std::cout << line;
            }
        }        

        if(regions != "-")
        {
            delete in;
        }

        if(out_bed != "" && out_bed != "-")
        {
            delete bed_out_stream;
        }

        if(out_errors != "" && out_errors != "-")
        {
            delete error_out_stream;
        }
        if(diff_out_stream)
        {
            delete diff_out_stream;
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

