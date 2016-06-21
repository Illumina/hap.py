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
 * \brief VCF Validator
 *
 * \details Validate VCF alleles and optionally expand haplotypes
 *
 * \file validatevcf.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "SimpleDiploidCompare.hh"
#include "Variant.hh"
#include "VariantInput.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"
#include "GraphReference.hh"
#include "DiploidCompare.hh"

#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>
#include <memory>

// error needs to come after boost headers.
#include "Error.hh"

/* #define DEBUG_VALIDATEVCF */

using namespace variant;
using namespace haplotypes;

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string ref_fasta;

    std::string chr = "";
    int64_t start = -1;
    int64_t end = -1;

    std::string file;
    std::string sample;

    std::string regions_bed = "";
    std::string targets_bed = "";

    std::string out_errors = "";

    // = max 12 unphased hets in segment
    int64_t blimit = -1;
    bool progress = false;
    int progress_seconds = 10;
    int max_n_haplotypes = 4096;
    int64_t hb_window = 30;
    int64_t hb_expand = 30;

    bool apply_filters = true;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-vcf", po::value<std::string>(), "VCF file to validate.")
            ("output-errors,e", po::value<std::string>(), "Output failure information in a bed file.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("location,l", po::value<std::string>(), "The location to start at.")
            ("regions,R", po::value<std::string>(), "Use a bed file for getting a subset of regions (traversal via tabix).")
            ("targets,T", po::value<std::string>(), "Use a bed file for getting a subset of targets (streaming the whole file, ignoring things outside the bed regions).")
            ("progress", po::value<bool>(), "Set to true to output progress information.")
            ("progress-seconds", po::value<int>(), "Output progress information every n seconds.")
            ("window,w", po::value<int64_t>(), "Overlap window to create haplotype blocks.")
            ("max-n-haplotypes,n", po::value<int>(), "Maximum number of haplotypes to enumerate.")
            ("expand-hapblocks", po::value<int64_t>(), "Number of bases to expand around each haplotype block.")
            ("limit", po::value<int64_t>(), "Maximum number of haplotype blocks to process.")
            ("apply-filters", po::value<bool>(), "Apply filtering in VCF (on by default).")
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
            std::cout << "validatevcf version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-vcf"))
        {
            std::string vr = vm["input-vcf"].as<std::string>();
            std::vector<std::string> v;
            stringutil::split(vr, v, ":");
            // in case someone passes a ":"
            assert(v.size() > 0);

            file = v[0];
            sample = "";
            if(v.size() > 1)
            {
                sample = v[1];
            }
        }

        if (vm.count("output-errors"))
        {
            out_errors = vm["output-errors"].as< std::string >();
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

        if (vm.count("expand-hapblocks"))
        {
            hb_expand = vm["expand-hapblocks"].as< int64_t >();
        }

        if (vm.count("window"))
        {
            hb_window = vm["window"].as< int64_t >();
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
            apply_filters = vm["apply-filters-truth"].as< bool >();
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
        vr.setReturnHomref(true);
        vr.setValidateRef(ref_fasta.c_str(), true);

        if(regions_bed != "")
        {
            vr.setRegions(regions_bed.c_str(), true);
        }
        if(targets_bed != "")
        {
            vr.setTargets(targets_bed.c_str(), true);
        }

        int r1 = vr.addSample(file.c_str(), sample.c_str());

        vr.setApplyFilters(apply_filters, r1);

        bool stop_after_chr_change = false;
        if(chr != "")
        {
            vr.rewind(chr.c_str(), start);
            stop_after_chr_change = true;
        }

        std::ostream * error_out_stream = NULL;
        if(out_errors == "-" || out_errors.empty())
        {
            out_errors = "-";
            error_out_stream = &std::cout;
        }
        else
        {
            error_out_stream = new std::ofstream(out_errors.c_str());
        }

        int64_t nhb = 0;
        int64_t last_pos = std::numeric_limits<int64_t>::max();

        // hap-block status + update
        std::list<Variants> block_variants;
        int64_t block_start = -1;
        int64_t block_end = -1;
        uint64_t sl_id = 0;

        GraphReference gr(ref_fasta.c_str());
        DiploidReference dr(gr);
        dr.setNPaths(max_n_haplotypes);

        uint64_t total_blocks = 0;
        uint64_t failed_blocks = 0;
        uint64_t import_failed_records = 0;

        const auto finish_block = [&block_variants, r1,
                                   &chr,
                                   &block_start,
                                   &block_end,
                                   &sl_id,
                                   &dr, max_n_haplotypes,
                                   hb_window,
                                   &error_out_stream,
                                   &total_blocks,
                                   &failed_blocks,
                                   &import_failed_records,
                                   hb_expand] () {
#ifdef DEBUG_VALIDATEVCF
            std::cerr << "FINISH block " << sl_id << " " << chr << ":" << block_start << "-" << block_end << "\n";
#endif
            sl_id++;
            std::list<std::string> result;
            int nhaps = 0;
            if(max_n_haplotypes > 0 && block_variants.size() > 0)
            {
#ifdef DEBUG_VALIDATEVCF
                for (auto vars : block_variants) {
                    std::cerr << "PROCESSED: " << vars << "\n";
                }
#endif
                try
                {
                    dr.setRegion(chr.c_str(), std::max((int64_t )0, block_start - hb_expand), block_end + hb_expand,
                                 block_variants, r1);
                    auto l = dr.result();
                    if(l.size() == 0)
                    {
                        result.push_back("DIPENUM_FAIL: unknown");
                    }
                    nhaps = (int)l.size();
                }
                catch (std::runtime_error const & e)
                {
                    nhaps = 0;
                    result.push_back(std::string("DIPENUM_FAIL: ") + e.what());
                }
            }

            bool block_is_failed = false;
            for(auto & vars : block_variants)
            {
                vars.calls[0] = vars.calls[r1];
                vars.calls.resize(1);

#ifdef DEBUG_VALIDATEVCF
                std::cerr << "RAW: " << vars << "\n";
#endif
                if((!vars.calls[0].isHomref() && !vars.calls[0].isNocall()) || vars.getInfoFlag("IMPORT_FAIL"))
                {
                    if(vars.getInfoFlag("IMPORT_FAIL"))
                    {
                        block_is_failed = true;
                        import_failed_records++;
                        result.push_back("IMPORT_FAIL");
                    }
                    if(result.size() > 0)
                    {
                        block_is_failed = true;
                        bool b = false;
                        std::string allerrors;
                        for (auto i : result) {
                            if(!allerrors.empty())
                            {
                                allerrors += ",";
                            }
                            allerrors += i;
                            b = true;
                        }
                        if(b)
                        {
                            allerrors += ",SUPERLOCUS_ID=" + std::to_string(sl_id);
                            allerrors += ",SUPERLOCUS_NHAPS=" + std::to_string(nhaps);
                            if(error_out_stream)
                            {
                                *error_out_stream << chr << "\t" << vars.pos << "\t"
                                                  << vars.pos + vars.len << "\t" << vars << "\t" << allerrors << "\n";
                            }
                        }
                    }
                }
            }
            if(block_is_failed)
            {
                ++failed_blocks;
            }
            ++total_blocks;
            block_variants.clear();
            block_start = -1;
            block_end = -1;
        };

        auto start_time = std::chrono::high_resolution_clock::now();
        auto last_time = std::chrono::high_resolution_clock::now();
        while(vr.advance())
        {
            if(blimit > 0 && nhb++ > blimit)
            {
                // reached record limit
                break;
            }
            Variants & v = vr.current();

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

            if (v.chr != chr || (block_end > 0 && block_end + hb_window < v.pos))
            {
                finish_block();
            }
            chr = v.chr;

            if (block_start < 0)
            {
                block_start = v.pos;
            }
            else
            {
                block_start = std::min(block_start, v.pos);
            }

            if (block_end < 0)
            {
                block_end = v.pos + v.len - 1;
            }
            else
            {
                block_end = std::max(v.pos + v.len - 1, block_end);
            }

            block_variants.push_back(v);

#ifdef DEBUG_VALIDATEVCF
            std::cerr << v << "\n";
            std::cerr << "block_start : " << block_start << "\t"
                      << "block_end : " << block_end << "\t"
                      << "block_size : " << block_variants.size() << "\t"
                      << "\n";
#endif

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
#ifdef DEBUG_VALIDATEVCF
        std::cerr << "END\n";
        std::cerr << "block_start : " << block_start << "\t"
                  << "block_end : " << block_end << "\t"
                  << "block_size : " << block_variants.size() << "\t"
                  << "\n";
#endif
        finish_block();

        std::cerr << "Total superloci: " << total_blocks <<
                     " failed: " << failed_blocks <<
                     " failed records: " << import_failed_records << "\n";

        if(error_out_stream && out_errors != "-")
        {
            delete error_out_stream;
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

