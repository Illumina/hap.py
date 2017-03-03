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
 * Output bed file to classify regions in a gVCF file
 *
 * \file quantify.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "Version.hh"
#include "Variant.hh"

#include "helpers/StringUtil.hh"

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <memory>
#include <htslib/synced_bcf_reader.h>
#include <helpers/BCFHelpers.hh>
#include <htslib/vcf.h>

#include "Error.hh"


using namespace variant;


int main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::string file;
    std::string output;
    std::string reference;
    std::string target;

    try
    {
        try
        {
            // Declare the supported options.
            po::options_description desc("Allowed options");
            desc.add_options()
                ("help,h", "produce help message")
                ("version", "Show version")
                ("input-file", po::value<std::string>(), "The input file")
                ("output-file,o", po::value<std::string>(), "The output file name (BED Format).")
                ("reference,r", po::value<std::string>(), "Reference fasta file.")
                ("target-region,T", po::value<std::string>(), "Optional bed file with target regions")
            ;

            po::positional_options_description popts;
            popts.add("input-file", -1);

            po::options_description cmdline_options;
            cmdline_options
                .add(desc);

            po::variables_map vm;

            po::store(po::command_line_parser(argc, argv).
                options(cmdline_options).positional(popts).run(), vm);
            po::notify(vm);

            if (vm.count("version"))
            {
                std::cout << "gvcf2bed version " << HAPLOTYPES_VERSION << "\n";
                return 0;
            }

            if (vm.count("help"))
            {
                std::cout << desc << "\n";
                return 1;
            }

            if (vm.count("input-file"))
            {
                file = vm["input-file"].as<std::string>();
            }

            if (vm.count("output-file"))
            {
                output = vm["output-file"].as<std::string>();
            }

            if (vm.count("reference"))
            {
                reference = vm["reference"].as<std::string>();
            }

            if (vm.count("target-region"))
            {
                target = vm["target-region"].as<std::string>();
            }

            if (file.empty())
            {
                std::cerr << "Please specify one input file / sample.\n";
                return 1;
            }

            if (output.empty())
            {
                std::cerr << "Please specify an output file.\n";
                return 1;
            }
        }
        catch (po::error &e)
        {
            std::cerr << e.what() << "\n";
            return 1;
        }

        FastaFile ref_fasta(reference.c_str());

        bcf_srs_t *reader = bcf_sr_init();
        reader->collapse = COLLAPSE_NONE;
            reader->require_index = 0;
            reader->streaming = 1;

        if(!target.empty())
        {
            bcf_sr_set_targets(reader, target.c_str(), 1, 0);
        }

        if (!bcf_sr_add_reader(reader, file.c_str()))
        {
            error("Failed to open or file not indexed: %s\n", file.c_str());
        }

        bcf_hdr_t *hdr = reader->readers[0].header;

        std::ostream *output_stream;
        std::unique_ptr<std::ofstream> output_stream_file;
        if (output == "-")
        {
            output_stream = &std::cout;
        }
        else
        {
            output_stream_file.reset(new std::ofstream(output));
            output_stream = output_stream_file.get();
        }

        struct VCFInterval
        {
            std::string chr;
            int64_t start = -1, end = -1;
            std::string extra;

            void print(std::ostream & o)
            {
                if(!chr.empty() && end >= 0)
                {
                    o << chr << "\t"
                      << start << "\t"
                      << end + 1 << "\t"
                      << extra
                      << std::endl;
                    chr = "";
                    start = end = -1;
                    extra = "";
                }
            }
        };

        VCFInterval current_interval;

        int nl = 1;
        while (nl)
        {
            nl = bcf_sr_next_line(reader);
            if (nl <= 0)
            {
                break;
            }
            if (!bcf_sr_has_line(reader, 0))
            {
                continue;
            }
            bcf1_t *line = reader->readers[0].buffer[0];

            bcf_unpack(line, BCF_UN_ALL);

            std::set<std::string> filters;

            for (int j = 0; j < line->d.n_flt; ++j)
            {
                std::string filter = "PASS";
                int k = line->d.flt[j];

                if (k >= 0)
                {
                    filter = bcf_hdr_int2id(hdr, BCF_DT_ID, line->d.flt[j]);
                }
                if (filter != "PASS" && !filter.empty())
                {
                    filters.insert(filter);
                }
            }

            const std::string vchr = bcfhelpers::getChrom(hdr, line);
            int64_t refstart, refend;
            bcfhelpers::getLocation(hdr, line, refstart, refend);

            const std::string ref_allele = line->d.allele[0];

            // pure insertion lines must be fully contained in CONF to match
            bool is_pure_insertion = false;

            if (bcfhelpers::classifyAlleleString(ref_allele).first == bcfhelpers::AlleleType::NUC)
            {
                is_pure_insertion = line->n_allele > 1;
                int64_t updated_ref_start = std::numeric_limits<int64_t>::max();
                int64_t updated_ref_end = std::numeric_limits<int64_t>::min();
                bool nuc_alleles = false;

                for (int al = 1; al < line->n_allele; ++al)
                {
                    RefVar al_rv;
                    al_rv.start = line->pos;
                    al_rv.end = (int64_t) (line->pos + ref_allele.size() - 1);
                    auto ca = bcfhelpers::classifyAlleleString(line->d.allele[al]);
                    if (ca.first == bcfhelpers::AlleleType::MISSING)
                    {
                        ca.second = "";
                    }
                    else if (ca.first != bcfhelpers::AlleleType::NUC)
                    {
                        break;
                    }
                    nuc_alleles = true;
                    al_rv.alt = ca.second;

                    variant::trimRight(ref_fasta, vchr.c_str(), al_rv, false);
                    variant::trimLeft(ref_fasta, vchr.c_str(), al_rv, false);

                    if (al_rv.end >= al_rv.start)
                    {
                        updated_ref_start = std::min(updated_ref_start, al_rv.start);
                        updated_ref_end = std::max(updated_ref_end, al_rv.end);
                        is_pure_insertion = false;
                    }
                    else
                    {
                        // this is an insertion *before* start, it will have end < start
                        // insertions are captured by the reference bases before and after
                        updated_ref_start = std::min(updated_ref_start, al_rv.start - 1);
                        updated_ref_end = std::max(updated_ref_end, al_rv.start);
                    }
                }

                if (nuc_alleles)
                {
                    refstart = updated_ref_start;
                    refend = updated_ref_end;
                }
            }

            static const int homref = 1;
            static const int call = 2;
            static const int nocall = 4;

            int gt_types = 0;

            for (int j = 0; j < line->n_sample; ++j)
            {
                int gt[MAX_GT]{-1, -1};
                int ngt = 0;
                bool phased = false;
                bcfhelpers::getGT(hdr, line, j, gt, ngt, phased);

                bool has_call = false;
                for (int g = 0; g < ngt; ++g)
                {
                    if (gt[g] > 0)
                    {
                        gt_types |= call;
                        has_call = true;
                    }
                }
                if (ngt == 0)
                {
                    gt_types |= nocall;
                }
                else if (!has_call)
                {
                    gt_types |= homref;
                }
            }

            std::string extra;
            auto addExtra = [&extra](std::string const &x)
            {
                if (!extra.empty())
                {
                    extra += ",";
                }
                extra += x;
            };

            if (gt_types & homref)
            {
                addExtra("homref");
            }
            if (gt_types & call)
            {
                addExtra("call");
            }
            if (gt_types & nocall)
            {
                addExtra("nocall");
            }
            if (is_pure_insertion)
            {
                addExtra("insertion");
            }
            if (extra.empty())
            {
                extra = "unknown";
            }

            extra += ":";
            if(filters.empty())
            {
                extra += "PASS";
            }
            else
            {
                int count = 0;
                for(const auto & f : filters)
                {
                    if(count > 0)
                    {
                        extra += ",";
                    }
                    extra += f;
                }
            }

            if( current_interval.chr == vchr &&
                refstart <= current_interval.end &&
                refend >= current_interval.start)
            {
                current_interval.start = std::min(refstart, current_interval.start);
                current_interval.end = std::max(refstart, current_interval.end);
                if(current_interval.extra.find(extra) == std::string::npos)
                {
                    current_interval.extra += ";" + extra;
                }
            }
            else
            {
                current_interval.print(*output_stream);
                current_interval.chr = vchr;
                current_interval.start = refstart;
                current_interval.end = refend;
                current_interval.extra = extra;
            }
        }

        current_interval.print(*output_stream);

        bcf_sr_destroy(reader);
    }
    catch (std::runtime_error &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (std::logic_error &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}

