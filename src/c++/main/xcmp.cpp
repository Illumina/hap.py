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
 *  \brief GT-GT comparison, followed by haplotype comparison for mismatching blocks
 *
 *
 * \file xcmp.cpp
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

/* #define DEBUG_XCMP */

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
    std::string file2;
    std::string sample2;

    std::string regions_bed = "";
    std::string targets_bed = "";

    std::string out_vcf = "";
    std::string out_errors = "";

    // = max 12 unphased hets in segment
    int64_t blimit = -1;
    bool progress = false;
    int progress_seconds = 10;
    int max_n_haplotypes = 4096;
    int64_t hb_window = 30;
    int64_t hb_expand = 30;

    bool apply_filters_query = false;
    bool apply_filters_truth = true;
    bool preprocess = false;
    bool leftshift = false;
    bool always_hapcmp = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-vcfs", po::value<std::vector<std::string> >(), "Two VCF files to compare (use file:sample for a specific sample column).")
            ("output-vcf,o", po::value<std::string>(), "Output variant comparison results to VCF.")
            ("output-errors,e", po::value<std::string>(), "Output failure information.")
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
            ("apply-filters-truth", po::value<bool>(), "Apply filtering in truth VCF (on by default).")
            ("apply-filters-query,f", po::value<bool>(), "Apply filtering in query VCF (off by default).")
            ("preprocess-variants,V", po::value<bool>(), "Apply variant normalisations, trimming, realignment for complex variants (off by default).")
            ("leftshift", po::value<bool>(), "Left-shift indel alleles (off by default).")
            ("always-hapcmp", po::value<bool>(), "Always compare haplotype blocks (even if they match). Testing use only/slow.")
        ;

        po::positional_options_description popts;
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
            std::cout << "xcmp version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
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

        if (vm.count("apply-filters-truth"))
        {
            apply_filters_truth = vm["apply-filters-truth"].as< bool >();
        }

        if (vm.count("apply-filters-query"))
        {
            apply_filters_query = vm["apply-filters-query"].as< bool >();
        }

        if (vm.count("always-hapcmp"))
        {
            always_hapcmp = vm["always-hapcmp"].as< bool >();
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

        if(regions_bed != "")
        {
            vr.setRegions(regions_bed.c_str(), true);
        }
        if(targets_bed != "")
        {
            vr.setTargets(targets_bed.c_str(), true);
        }

        int r1 = vr.addSample(file1.c_str(), sample1.c_str());
        int r2 = vr.addSample(file2.c_str(), sample2.c_str());

        vr.setApplyFilters(apply_filters_truth, r1);
        vr.setApplyFilters(apply_filters_query, r2);

        VariantInput vi(
            ref_fasta.c_str(),
            preprocess || leftshift,          // bool leftshift
            preprocess || leftshift,          // bool refpadding
            true,                // bool trimalleles = false, (remove unused alleles)
            preprocess || leftshift,      // bool splitalleles = false,
            ( preprocess || leftshift ) ? 2 : 0,  // int mergebylocation = false,
            true,                // bool uniqalleles = false,
            true,                // bool calls_only = true,
            false,               // bool homref_split = false // this is handled by calls_only
            preprocess,          // bool primitives = false
            false,               // bool homref_output
            leftshift ? hb_window-1 : 0   // int64_t leftshift_limit
            );

        VariantProcessor & vp = vi.getProcessor();

        vp.setReader(vr, VariantBufferMode::buffer_block, 10*hb_window);

        bool stop_after_chr_change = false;
        if(chr != "")
        {
            vp.rewind(chr.c_str(), start);
            stop_after_chr_change = true;
        }

        std::unique_ptr<VariantWriter> pvw;
        if (out_vcf != "")
        {
            pvw = std::move(std::unique_ptr<VariantWriter> (new VariantWriter(out_vcf.c_str(), ref_fasta.c_str())));
            pvw->addHeader(vr);
            pvw->addHeader("##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">");
            pvw->addHeader("##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">");
            pvw->addHeader("##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
            pvw->addHeader("##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
            pvw->addHeader("##INFO=<ID=ctype,Number=1,Type=String,Description=\"Type of comparison performed\">");
            pvw->addHeader("##INFO=<ID=HapMatch,Number=0,Type=Flag,Description=\"Variant is in matching haplotype block\">");
            pvw->addSample("TRUTH");
            pvw->addSample("QUERY");
        }

        std::ostream * error_out_stream = NULL;
        if(out_errors == "-")
        {
            error_out_stream = &std::cerr;
        }
        else if(out_errors != "")
        {
            error_out_stream = new std::ofstream(out_errors.c_str());
        }
        DiploidCompare hc(ref_fasta.c_str());
        hc.setMaxHapEnum(max_n_haplotypes);
        hc.setDoAlignments(false);

        int64_t nhb = 0;
        int64_t last_pos = std::numeric_limits<int64_t>::max();

        // hap-block status + update
        std::list<Variants> block_variants;
        int64_t block_start = -1;
        int64_t block_end = -1;
        int n_nonsnp = 0, calls_1 = 0, calls_2 = 0;
        bool has_mismatch = false;

        const auto finish_block = [&block_variants, r1, r2,
                                   &chr,
                                   &block_start,
                                   &block_end,
                                   &n_nonsnp, &calls_1, &calls_2,
                                   &has_mismatch,
                                   &pvw, &error_out_stream,
                                   &hc,
                                   hb_expand, always_hapcmp] ()
        {
            bool hap_match = false, hap_fail = false, hap_run = false;
            // try HC if we have mismatches, and if the number of calls is > 0
            if (always_hapcmp || (has_mismatch && calls_1 > 0 && calls_2 > 0 && n_nonsnp > 0))
            {
                try
                {
                    hap_run = true;
                    hap_fail = true;
                    hc.setRegion(chr.c_str(), std::max(int64_t(0), block_start-hb_expand), block_end + hb_expand,
                                 block_variants, r1, r2);
                    DiploidComparisonResult const & hcr = hc.getResult();
#ifdef DEBUG_XCMP
                    std::cerr << hcr << "\n";
#endif
                    hap_match = hcr.outcome == dco_match;
                    hap_fail = !(hcr.outcome == dco_match || hcr.outcome == dco_mismatch);
                }
                catch(std::runtime_error &e)
                {
                    if (error_out_stream)
                    {
                        *error_out_stream << chr << "\t" << block_start << "\t" << block_end+1 << "\t" << "hap_error\t" << e.what() << "\n";
                    }
                }
                catch(std::logic_error &e)
                {
                    if (error_out_stream)
                    {
                        *error_out_stream << chr << "\t" << block_start << "\t" << block_end+1 << "\t" << "hap_error\t" << e.what() << "\n";
                    }
                }
            }
            std::string result;
            if(hap_run)
            {
                if (hap_fail)
                {
                    result = "hapfail:";
                }
                else if(has_mismatch)
                {
                    result = "hap:";
                }
                else
                {
                    result = "simple:";
                }
            }
            else
            {
                result = "simple:";
            }

            if(hap_run && !has_mismatch && !hap_match)
            {
                result += "suspicious_simple_match";
            }
            else if(always_hapcmp && hap_match && ((calls_1 == 0 && calls_2 > 0) || (calls_1 > 0 && calls_2 == 0)))
            {
                bool any_filtered = false;
                for (Variants const & v : block_variants)
                {
                    for (Call const & c : v.calls)
                    {
                        for (size_t i = 0; i < c.nfilter; ++i)
                        {
                            if(c.filter[i] != "PASS" && c.filter[i] != ".")
                            {
                                any_filtered = true;
                                break;
                            }
                        }
                    }
                    if(any_filtered)
                    {
                        break;
                    }
                }

                if(any_filtered)
                {
                    result += "match_ignoring_filtered";
                }
                else
                {
                    result += "suspicious_hap_match";
                }
            }
            else if(hap_match || !has_mismatch)
            {
                result += "match";
            }
            else
            {
                result += "mismatch";
            }

            if(error_out_stream)
            {
                *error_out_stream << chr << "\t" << block_start << "\t" << block_end+1 << "\t" << result << "\t"
                                  << has_mismatch << ":" << hap_match << ":" << hap_fail << ":"
                                  << calls_1 << ":" << calls_2 << ":" << n_nonsnp << "\n";
            }
            if (pvw)
            {
                for (Variants & v : block_variants)
                {
                    if(v.info != "")
                    {
                        v.info += ";";
                    }
                    v.info += std::string("ctype=") + result;
                    if (hap_match)
                    {
                        v.info += ";HapMatch";
                    }
                    pvw->put(v);
                }
            }

            block_variants.clear();
            block_start = -1;
            block_end = -1;
            n_nonsnp = 0;
            calls_1 = 0;
            calls_2 = 0;
            has_mismatch = false;
        };

        auto start_time = std::chrono::high_resolution_clock::now();
        auto last_time = std::chrono::high_resolution_clock::now();
        while(vp.advance())
        {
            if(blimit > 0 && nhb++ > blimit)
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

            if(compareVariants(v, r1, r2, n_nonsnp, calls_1, calls_2) != dco_match)
            {
                has_mismatch = true;
            }
            block_variants.push_back(v);

#ifdef DEBUG_XCMP
            std::cerr << v << "\n";
            std::cerr << "block_start : " << block_start << "\t"
                      << "block_end : " << block_end << "\t"
                      << "block_size : " << block_variants.size() << "\t"
                      << "n_nonsnp : " << n_nonsnp << "\t"
                      << "calls_1 : " << calls_1 << "\t"
                      << "calls_2 : " << calls_2 << "\t"
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
#ifdef DEBUG_XCMP
        std::cerr << "END\n";
        std::cerr << "block_start : " << block_start << "\t"
                  << "block_end : " << block_end << "\t"
                  << "block_size : " << block_variants.size() << "\t"
                  << "n_nonsnp : " << n_nonsnp << "\t"
                  << "calls_1 : " << calls_1 << "\t"
                  << "calls_2 : " << calls_2 << "\t"
                  << "\n";
#endif
        finish_block();
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

