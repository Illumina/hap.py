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
#include <queue>
#include <mutex>
#include <future>
#include <htslib/synced_bcf_reader.h>
#include <helpers/BCFHelpers.hh>
#include <htslib/vcf.h>

#include "Error.hh"

#include "BlockQuantify.hh"

using namespace variant;


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::string file;
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

    int threads = 1;
    int blocksize = 20000;

    QuantifyRegions regions;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::string>(), "The input file")
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
            ("threads", po::value<int>(), "Number of threads to use.")
            ("blocksize", po::value<int>(), "Number of variants per block.")
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
            file = vm["input-file"].as< std::string >();
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

        if (vm.count("threads"))
        {
            threads = vm["threads"].as< int >();
        }

        if (vm.count("blocksize"))
        {
            blocksize = vm["blocksize"].as< int >();
        }

        if(file.size() == 0)
        {
            std::cerr << "Please specify one input file / sample.\n";
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
        bcf_srs_t * reader = bcf_sr_init();
        if(!only_regions.empty()) {
            int result = bcf_sr_set_regions(reader, only_regions.c_str(), 1);
            if(result < 0)
            {
                error("Failed to set regions string %s.", only_regions.c_str());
            }
        }
        if (!bcf_sr_add_reader(reader, file.c_str()))
        {
            error("Failed to open or file not indexed: %s\n", file.c_str());
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
            }
            if(success < 0)
            {
                error("Cannot seek to %s:%i", chr.c_str(), start);
            }
        }

        bcf_hdr_t * hdr = reader->readers[0].header;
        htsFile * writer = nullptr;

        if (output_vcf != "")
        {
            bcf_hdr_append(hdr, "##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">");
            bcf_hdr_append(hdr, "##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">");
            bcf_hdr_append(hdr, "##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
            bcf_hdr_append(hdr, "##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
            bcf_hdr_append(hdr, "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">");
            bcf_hdr_append(hdr, "##INFO=<ID=T_VT,Number=1,Type=String,Description=\"High-level variant type in truth (SNP|INDEL).\">");
            bcf_hdr_append(hdr, "##INFO=<ID=Q_VT,Number=1,Type=String,Description=\"High-level variant type in query (SNP|INDEL).\">");
            bcf_hdr_append(hdr, "##INFO=<ID=T_LT,Number=1,Type=String,Description=\"High-level location type in truth (het|hom|hetalt).\">");
            bcf_hdr_append(hdr, "##INFO=<ID=Q_LT,Number=1,Type=String,Description=\"High-level location type in query (het|hom|hetalt).\">");
            if(output_vtc)
            {
                bcf_hdr_append(hdr, "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">");
            }

            const char * mode = "wu";

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
            bcf_hdr_write(writer, hdr);
        }

        /** local function to count variants in all samples */
        int64_t rcount = 0;
        std::string current_chr = "";
        int vars_in_block = 0;
        std::unique_ptr<BlockQuantify> p_bq(new BlockQuantify(hdr, regions, ref_fasta, output_vtc, count_homref));

        /** async stuff. each block can be counted in parallel, but we need to
         *  write out the variants sequentially.
         *  Therefore, we keep a future for each block to be able to join
         *  when it's processed
         */
        std::queue<std::pair <
            std::future<void>,
            std::unique_ptr<BlockQuantify>
        >> blocks;

        std::map<std::string, VariantStatistics> count_map;
        auto output_counts = [&count_map, &writer, &blocks, hdr](int min_size) {
            while(blocks.size() > (unsigned )min_size)
            {
                blocks.front().first.get();
                if(writer)
                {
                    auto const & variants = blocks.front().second->getVariants();
                    for(auto & v : variants)
                    {
                        bcf_write1(writer, hdr, v);
                    }
                }

                auto const & cm = blocks.front().second->getCounts();
                for(auto const & c : cm)
                {
                    auto it = count_map.find(c.first);
                    if (it == count_map.end()) {
                        count_map.insert(c);
                    }
                    else
                    {
                        it->second.add(c.second);
                    }
                }

                blocks.pop();
            }
        };

        int nl = 1;
        while(nl)
        {
            nl = bcf_sr_next_line(reader);
            if (nl <= 0)
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
            const std::string vchr = bcfhelpers::getChrom(hdr, line);
            if(end != -1 && ((!current_chr.empty() && vchr != current_chr) || line->pos > end))
            {
                break;
            }

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
                        filter = bcf_hdr_int2id(hdr, BCF_DT_ID, line->d.flt[j]);
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

            current_chr = vchr;

            p_bq->add(bcf_dup(line));

            ++vars_in_block;
            if(vars_in_block > blocksize)
            {
                std::future<void> f = std::async(std::launch::async, &BlockQuantify::count, p_bq.get());
                // clear / write out some blocks (make sure we have at least 2xthreads tasks left)
                output_counts(threads);
                blocks.emplace(std::move(f), std::move(p_bq));
                p_bq.reset(new BlockQuantify(hdr, regions, ref_fasta, output_vtc, count_homref));
                vars_in_block = 0;
            }

            if (message > 0 && (rcount % message) == 0)
            {
                std::cout << stringutil::formatPos(vchr.c_str(), line->pos) << "\n";
            }
            // count variants here
            ++rcount;
        }

        {
            std::future<void> f = std::async(&BlockQuantify::count, p_bq.get());
            // clear / write out some blocks (make sure we have at least 2xthreads tasks left)
            blocks.emplace(std::move(f), std::move(p_bq));
        }
        // clear remaining
        output_counts(0);

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

        hts_close(writer);
        bcf_sr_destroy(reader);
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
