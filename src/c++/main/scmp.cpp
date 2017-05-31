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
 * Somatic / by allele comparison -- checks if we see the same alleles (ignoring genotypes).
 *
 * \file scmp.cpp
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
#include <queue>
#include <mutex>
#include <future>
#include <htslib/synced_bcf_reader.h>
#include <helpers/BCFHelpers.hh>
#include <htslib/vcf.h>

#include "Error.hh"
#include "helpers/Timing.hh"

#include "BlockAlleleCompare.hh"

using namespace variant;

/* #define DEBUG_SCMP_THREADING */

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    try
    {
        std::string input_vcf;
        std::string output_vcf;
        std::string ref;
        std::string only_regions;
        std::string qq_field = "QUAL";
        std::string comparison_mode;

        // limits
        std::string chr;
        int64_t start = -1;
        int64_t end = -1;
        int64_t rlimit = -1;
        int64_t message = -1;
        bool apply_filters = false;
        int threads = 1;
        int blocksize = 20000;
        int min_var_distance = 100000;

        // comparison mode
        BlockAlleleCompare::ComparisonMode cmode = BlockAlleleCompare::ALLELES;
        Json::Value cmode_params;

        try
        {
            // Declare the supported options.
            po::options_description desc("Allowed options");
            desc.add_options()
                    ("help,h", "produce help message")
                    ("version", "Show version")
                    ("input-file", po::value<std::string>(), "Input VCF file. Must have exactly two samples, the first "
                            "sample will be used as truth, the second one as query. This can be obtained using bcftools: "
                            "bcftools merge truth.vcf.gz query.vcf.gz --force-samples")
                    ("output-file,o", po::value<std::string>(), "The output file name (VCF / BCF / VCF.gz).")
                    ("reference,r", po::value<std::string>(), "The reference fasta file (needed only for VCF output).")
                    ("location,l", po::value<std::string>(), "Start location.")
                    ("qq-field,qq", po::value<std::string>(),
                     "QQ field name -- this can be QUAL, an INFO or FORMAT field name")
                    ("only,O", po::value<std::string>(), "Bed file of locations (equivalent to -R in bcftools)")
                    ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
                    ("message-every", po::value<int64_t>(), "Print a message every N records.")
                    ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
                    ("threads", po::value<int>(), "Number of threads to use.")
                    ("blocksize", po::value<int>(), "Number of variants per block.")
                    ("min-var-distance", po::value<int>(), "Minimum distance between variants allowed to end up "
                            "in separate blocks")
                    ("comparison-mode,M", po::value<std::string>()->default_value("allele"),
                     "How to compare variants: allele (default) / distance / enumerate")
                    ("distance-maxdist", po::value<int>()->default_value(50),
                     "For distance comparison, this is the maximum distance between variants s.t. they get matched.")
                    ("enumerate-maxenum", po::value<int>()->default_value(1 << 16),
                     "For haplotype enumeration comparison, this is the number of possibilities to enumerate before giving up")
                    ;

            po::positional_options_description popts;
            popts.add("input-file", 1);

            po::options_description cmdline_options;
            cmdline_options.add(desc);

            po::variables_map vm;

            po::store(po::command_line_parser(argc, argv).
                    options(cmdline_options).positional(popts).run(), vm);
            po::notify(vm);

            if(vm.count("version"))
            {
                std::cout << "scmp version " << HAPLOTYPES_VERSION << "\n";
                return 0;
            }

            if(vm.count("help"))
            {
                std::cout << desc << "\n";
                return 1;
            }

            if(vm.count("input-file"))
            {
                input_vcf = vm["input-file"].as<std::string>();
            }
            else
            {
                error("Input file is required.");
            }

            if(vm.count("output-file"))
            {
                output_vcf = vm["output-file"].as<std::string>();
            }
            else
            {
                error("Output file name is required.");
            }

            if(vm.count("reference"))
            {
                ref = vm["reference"].as<std::string>();
            }
            else
            {
                error("To write an output VCF, you need to specify a reference file, too.");
            }

            if(vm.count("location"))
            {
                stringutil::parsePos(vm["location"].as<std::string>(), chr, start, end);
            }

            if(vm.count("qq-field"))
            {
                qq_field = vm["qq-field"].as<std::string>();
            }

            if(vm.count("only"))
            {
                only_regions = vm["only"].as<std::string>();
            }

            if(vm.count("limit-records"))
            {
                rlimit = vm["limit-records"].as<int64_t>();
            }

            if(vm.count("message-every"))
            {
                message = vm["message-every"].as<int64_t>();
            }

            if(vm.count("apply-filters"))
            {
                apply_filters = vm["apply-filters"].as<bool>();
            }

            if(vm.count("threads"))
            {
                threads = vm["threads"].as<int>();
            }

            if(vm.count("blocksize"))
            {
                blocksize = vm["blocksize"].as<int>();
            }

            if(vm.count("min-var-distance"))
            {
                min_var_distance = vm["min-var-distance"].as<int>();
            }

            comparison_mode = vm["comparison-mode"].as<std::string>();
            // make sure we don't break distance-based matching by blocking
            if(comparison_mode == "distance")
            {
                const int distance_max_dist = vm["distance-maxdist"].as<int>();
                min_var_distance = std::max(min_var_distance, distance_max_dist + 2);
                // ensure minimum block size also
                blocksize = std::max(blocksize, 20000);

                cmode = BlockAlleleCompare::DISTANCE;
                cmode_params["max_distance"] = distance_max_dist;
            }
            else if(comparison_mode == "alleles")
            {
                cmode = BlockAlleleCompare::ALLELES;
            }
            else if(comparison_mode == "enumerate")
            {
                cmode = BlockAlleleCompare::ENUMERATE_DIPLOID;
                const int max_enum = vm["distance-maxdist"].as<int>();
                cmode_params["max_enum"] = max_enum;
            }
            else
            {
                error("Invalid comparison mode: %s", comparison_mode.c_str());
            }
        }
        catch(po::error &e)
        {
            std::cerr << e.what() << "\n";
            return 1;
        }

        FastaFile ref_fasta(ref.c_str());
        bcf_srs_t *reader = bcf_sr_init();
        reader->collapse = COLLAPSE_NONE;

        if(!chr.empty() || !only_regions.empty())
        {
            reader->require_index = 1;
            reader->streaming = 0;
        }
        else
        {
            reader->require_index = 0;
            reader->streaming = 1;
        }

        if(!only_regions.empty())
        {
            int result = bcf_sr_set_regions(reader, only_regions.c_str(), 1);
            if(result < 0)
            {
                error("Failed to set regions string %s.", only_regions.c_str());
            }
        }
        if(!bcf_sr_add_reader(reader, input_vcf.c_str()))
        {
            error("Failed to open or file not indexed: %s\n", input_vcf.c_str());
        }

        if(!chr.empty())
        {
            int success = 0;
            if(start < 0)
            {
                success = bcf_sr_seek(reader, chr.c_str(), 0);
            }
            else
            {
                success = bcf_sr_seek(reader, chr.c_str(), start);
#ifdef DEBUG_SCMP
                std::cerr << "starting at " << chr << ":" << start << "\n";
#endif
            }
            if(success <= -2)
            {
                error("Cannot seek to %s:%i", chr.c_str(), start);
            }
        }

        auto hdr = bcfhelpers::ph(bcf_hdr_dup(reader->readers[0].header));

        BlockAlleleCompare::updateHeader(hdr.get());

        // rename the samples in the output header
        bcfhelpers::p_bcf_hdr output_header(bcfhelpers::ph(bcf_hdr_init("w")));
        {
            int len = 0;
            char *hdr_text = bcf_hdr_fmt_text(hdr.get(), 0, &len);
            if(!hdr_text)
            {
                error("Failed to process input VCF header.");
            }

            std::vector<std::string> split_header;
            stringutil::split(std::string(hdr_text, (unsigned long) len), split_header, "\n");
            free(hdr_text);
            for(std::string hl : split_header)
            {
                bcf_hdr_append(output_header.get(), hl.c_str());
            }
            bcf_hdr_add_sample(output_header.get(), "TRUTH");
            bcf_hdr_add_sample(output_header.get(), "QUERY");
            bcf_hdr_sync(output_header.get());
        }

        htsFile *writer = nullptr;

        const char *mode = "wu";

        if(stringutil::endsWith(output_vcf, ".vcf.gz"))
        {
            mode = "wz";
        }
        else if(stringutil::endsWith(output_vcf, ".bcf"))
        {
            mode = "wb";
        }

        if(!output_vcf.empty() && output_vcf[0] == '-')
        {
            writer = hts_open("-", mode);
        }
        else
        {
            writer = hts_open(output_vcf.c_str(), mode);
        }
        bcf_hdr_write(writer, output_header.get());

        /** local function to count variants in all samples */
        int64_t rcount = 0;
        std::string current_chr = "";
        int vars_in_block = 0;

        /*
         * we keep a list of jobs as futures which get processed by a pool of
         * workers
         */
        std::queue<std::unique_ptr<BlockAlleleCompare> > blocks;
        std::mutex blocks_mutex;

        /*
         * once a worker completes a comparison, the result goes into the list of
         * processed blocks.
         */
        std::list<std::unique_ptr<BlockAlleleCompare>> processed_blocks;
        std::mutex processed_blocks_mutex;


        /**
         * gets set to true by main thread
         */
        bool all_blocks_added = false;

        std::vector<std::thread> workers;

        auto worker = [&blocks,
                &blocks_mutex,
                &processed_blocks,
                &processed_blocks_mutex,
                &all_blocks_added]
        {
            bool work_left = true;
            while(work_left)
            {
                std::unique_ptr<BlockAlleleCompare> current_block;
#ifdef DEBUG_SCMP_THREADING
                std::cerr << "Worker " << std::thread::id() << " waiting." << std::endl;
#endif
                {
                    std::lock_guard<std::mutex> lk(blocks_mutex);
                    /* std::unique_lock<std::mutex> lk(blocks_mutex); */
                    /* lk.wait([&blocks]() { !blocks.empty(); }); */
                    if(!blocks.empty())
                    {
                        current_block = std::move(blocks.front());
                        blocks.pop();
                    }
                }
                if(current_block)
                {
#ifdef DEBUG_SCMP_THREADING
                    std::cerr << "Worker " << std::thread::id() << " has a job." << std::endl;
#endif
                    current_block->run();
#ifdef DEBUG_SCMP_THREADING
                    std::cerr << "Worker " << std::thread::id() << " finished a job." << std::endl;
#endif
                    {
                        std::lock_guard<std::mutex> lock(processed_blocks_mutex);
                        processed_blocks.emplace_back(std::move(current_block));
#ifdef DEBUG_SCMP_THREADING
                        std::cerr << "Worker " << std::thread::id() << " done." << std::endl;
#endif
                    }
                }
                work_left = !all_blocks_added || !blocks.empty();
            }
        };

        for(int j = 0; j < std::max(1, threads); ++j)
        {
            workers.emplace_back(worker);
        }

        std::unique_ptr<BlockAlleleCompare> p_bac(new BlockAlleleCompare(hdr, ref_fasta, qq_field));
        p_bac->setComparisonMode(cmode);
        p_bac->setComparisonParameters(cmode_params);

        int nl = 1;
        int previous_pos = -1;
        auto t0 = CPUClock::now();
        while(nl)
        {
            nl = bcf_sr_next_line(reader);
            if(nl <= 0)
            {
                break;
            }
            if(!bcf_sr_has_line(reader, 0))
            {
                continue;
            }
            bcf1_t *line = reader->readers[0].buffer[0];

            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }
            const std::string vchr = bcfhelpers::getChrom(hdr.get(), line);
            if(end != -1 && ((!current_chr.empty() && vchr != current_chr) || line->pos > end))
            {
                break;
            }

            if(!current_chr.empty() && vchr != current_chr)
            {
                // reset bs on chr switch
                previous_pos = -1;
                if(vars_in_block > 0)
                {
                    // make sure we create a new block
                    vars_in_block = blocksize + 1;
                }
            }
            current_chr = vchr;

            if(apply_filters)
            {
                bcf_unpack(line, BCF_UN_FLT);

                bool fail = false;
                for(int j = 0; j < line->d.n_flt; ++j)
                {
                    std::string filter = "PASS";
                    int k = line->d.flt[j];

                    if(k >= 0)
                    {
                        filter = bcf_hdr_int2id(hdr.get(), BCF_DT_ID, line->d.flt[j]);
                    }
                    if(filter != "PASS")
                    {
                        fail = true;
                        break;
                    }
                }

                // skip failing
                if(fail)
                {
                    continue;
                }
            }

            const int dist_to_previous = (previous_pos < 0) ? 0 : line->pos - previous_pos;
            previous_pos = line->pos;

            if(vars_in_block > blocksize || dist_to_previous > min_var_distance)
            {
                {
                    std::lock_guard<std::mutex> lock(blocks_mutex);
                    blocks.emplace(std::move(p_bac));
                }
                p_bac = std::move(
                        std::unique_ptr<BlockAlleleCompare>(new BlockAlleleCompare(hdr, ref_fasta, qq_field)));
                p_bac->setComparisonMode(cmode);
                p_bac->setComparisonParameters(cmode_params);
                vars_in_block = 0;
                previous_pos = -1;
            }

            p_bac->add(line);
            ++vars_in_block;

            if(message > 0 && (rcount % message) == 0)
            {
                auto t1 = CPUClock::now();
                typedef std::chrono::duration<double, typename CPUClock::period> Cycle;
                std::cout << stringutil::formatPos(vchr.c_str(), line->pos) << " cycles/record: " <<
                          (Cycle(t1 - t0).count() / rcount) << "\n";
            }
            // count variants here
            ++rcount;
        }

        {
            std::lock_guard<std::mutex> lock(blocks_mutex);
            blocks.emplace(std::move(p_bac));
            all_blocks_added = true;
        }

        if(message > 0)
        {
            std::cout << "All records read. Waiting for workers to complete" << std::endl;
        }

        for(auto &t : workers)
        {
            t.join();
        }

        // make sure we write in the correct order
        processed_blocks.sort([](std::unique_ptr<BlockAlleleCompare> const &b1,
                                 std::unique_ptr<BlockAlleleCompare> const &b2)
                              {
                                  return *b1 < *b2;
                              });

        for(auto &result : processed_blocks)
        {
            result->output(writer);
        }

        hts_close(writer);
        bcf_sr_destroy(reader);
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
