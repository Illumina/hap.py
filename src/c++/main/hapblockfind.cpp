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
 * \brief Find VCF regions that have overlapping variation within a fixed window
 *
 * \file hapblockfind.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "Variant.hh"
#include "Haplotype.hh"
#include "GraphReference.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"

#include <fstream>

// error needs to come after program_options. 
#include "Error.hh"

using namespace variant;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string ref_fasta;
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
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("limit-records,L", po::value<int64_t>(), "Maximum number of records to process")
            ("message-every,m", po::value<int64_t>(), "Print a message every N records.")
            ("window,w", po::value<int64_t>(), "Overlap window length.")
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
            std::cout << "hapblockfind version " << HAPLOTYPES_VERSION << "\n";
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
        else
        {
            error("Please specify a position / chromosome name using -l.");
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

        if(regions_bed != "")
        {
            r.setRegions(regions_bed.c_str(), true);
        }
        if(targets_bed != "")
        {
            r.setTargets(targets_bed.c_str(), true);
        }

        std::list<int> sids;
        for(size_t i = 0; i < files.size(); ++i)
        {
            sids.push_back(r.addSample(files[i].c_str(), samples[i].c_str()));
        }

        r.setApplyFilters(apply_filters);

        if(chr != "")
        {
            r.rewind(chr.c_str(), start);
        }

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

        std::string block_chr = "";
        int64_t rcount = 0;
        int64_t last_end = -1;

        int64_t block_start = -1;
        int hets = 0;
        int vars = 0;
        std::vector<size_t> present;
        present.resize(sids.size(), 0);
        std::vector<size_t> unphased_het;
        unphased_het.resize(sids.size(), 0);

        bool cont = true;
        while(cont && r.advance())
        {
            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    cont = false;
                    break;
                }
            }

            Variants & v = r.current();
            if((end != -1 && v.pos > end) || (chr != "" && v.chr != chr))
            {
                cont = false;
                break;
            }
            chr = v.chr;

            if(message > 0 && (rcount % message) == 0)
            {
                std::cerr << v << "\n";
            }

            int hets_this_pos = 0;
            size_t homref = 0;

            for(int s : sids)
            {
                Call & c = v.calls[s];
                gttype gtt = getGTType(c);
                // count unphased hets only
                if((gtt == gt_het || gtt == gt_hetalt) && !c.phased)
                {
                    unphased_het[s]++;
                    ++hets_this_pos;
                }
                if(gtt == gt_homref || gtt == gt_unknown)
                {
                    ++homref;
                }
            }

            if(homref == sids.size())
            {
                continue;
            }

            if((last_end < 0 || v.pos <= last_end + window) && cont)
            {
                block_chr = v.chr;
                // extend block
                hets += hets_this_pos;
                vars++;
                last_end = std::max(last_end, v.pos + v.len - 1);

                if (block_start < 0)
                {
                    block_start = v.pos;
                }
            }
            else if(cont)
            {
                if(vars > 0)
                {
                    // finish block
                    (*outputfile) << block_chr << "\t" << std::max((int64_t)0, block_start) 
                               << "\t" << last_end + 1
                               << "\t" << vars << "\t" << hets
                               << "\n";                    
                    if(verbose)
                    {
                        std::cerr << "Finished block of length " << (1 + last_end - block_start) << " at " << v.chr << ":" << v.pos << "\n";
                    }
                }

                block_chr = v.chr;
                block_start = v.pos;
                hets = hets_this_pos;
                vars = 1;
                last_end = v.pos + v.len - 1;
                if(verbose)
                {
                    std::cerr << "Starting new block after " << v.chr << ":"  << block_start << "-" << last_end << "\n";
                }
            }

            ++rcount;
        }

        if(vars > 0)
        {
            // finish block
            (*outputfile) << block_chr << "\t" << std::max((int64_t)0, block_start) 
                       << "\t" << last_end + 1
                       << "\t" << vars << "\t" << hets
                       << "\n";                    
            if(verbose)
            {
                std::cerr << "[END] Finished block of length " << (1 + last_end - block_start) << " at " << block_chr << ":" << last_end << "\n";
            }
        }

        if(out_bed != "-" && out_bed != "")
        {
            delete outputfile;
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

