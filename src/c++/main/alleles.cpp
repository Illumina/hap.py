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
 * Somatic VCF preprocessor
 *
 * \file alleles.cpp
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

using namespace variant;


int main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    try
    {
        std::string input_vcf;
        std::string output_vcf;
        std::string sample = "SAMPLE";
        enum
        {
            HEMI, HET, HOM, HALF, FIRST
        } gt_output = HALF;
        try
        {
            // Declare the supported options.
            po::options_description desc("Allowed options");
            desc.add_options()
                    ("help,h", "produce help message")
                    ("version", "Show version")
                    ("input-file", po::value<std::string>(), "Input VCF file.")
                    ("output-file,o", po::value<std::string>(), "The output file name (VCF / BCF / VCF.gz).")
                    ("gt", po::value<std::string>(),
                     "What GT to write: hemi | het | hom | half | first for 1 | 0/1 | 1/1 | ./1 | "
                             "first sample GT; "
                             "default is half")
                    ("sample", po::value<std::string>()->default_value("SAMPLE"), "Output sample name");

            po::positional_options_description popts;
            popts.add("input-file", 1);

            po::options_description cmdline_options;
            cmdline_options
                    .add(desc);

            po::variables_map vm;

            po::store(po::command_line_parser(argc, argv).
                    options(cmdline_options).positional(popts).run(), vm);
            po::notify(vm);

            if(vm.count("version"))
            {
                std::cout << "alleles version " << HAPLOTYPES_VERSION << "\n";
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

            if(vm.count("gt"))
            {
                const std::string gt_str = vm["gt"].as<std::string>();

                if(gt_str == "het")
                {
                    gt_output = HET;
                }
                else if(gt_str == "hom")
                {
                    gt_output = HOM;
                }
                else if(gt_str == "hemi")
                {
                    gt_output = HEMI;
                }
                else if(gt_str == "half")
                {
                    gt_output = HALF;
                }
                else if(gt_str == "first")
                {
                    gt_output = FIRST;
                }
                else
                {
                    error("Unknown GT output type: %s", gt_str.c_str());
                }
            }

            sample = vm["sample"].as<std::string>();
        }
        catch(po::error &e)
        {
            std::cerr << e.what() << "\n";
            return 1;
        }

        bcf_srs_t *reader = bcf_sr_init();
        reader->collapse = COLLAPSE_NONE;
        reader->require_index = 0;
        reader->streaming = 1;

        if(!bcf_sr_add_reader(reader, input_vcf.c_str()))
        {
            error("Failed to open or file not indexed: %s\n", input_vcf.c_str());
        }

        bcf_hdr_t *hdr = reader->readers[0].header;

        if(!bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "GT", NULL))
        {
            bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        }

        bcf_hdr_sync(hdr);

        // rename the samples in the output header
        bcfhelpers::p_bcf_hdr output_header(bcfhelpers::ph(bcf_hdr_init("w")));
        {
            int len = 0;
            char *hdr_text = bcf_hdr_fmt_text(hdr, 0, &len);
            if(!hdr_text)
            {
                error("Failed to process input VCF header.");
            }
            std::vector<std::string> split_header;
            stringutil::split(std::string(hdr_text, (unsigned long) len), split_header, "\n");
            free(hdr_text);
            for(std::string hl : split_header)
            {
                // formats are appended below
                bcf_hdr_append(output_header.get(), hl.c_str());
            }
            bcf_hdr_append(output_header.get(), "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            bcf_hdr_add_sample(output_header.get(), sample.c_str());
        }

        std::set<std::string> format_names;
        for(int i = 0; i < hdr->nhrec; i++)
        {
            const auto hr = hdr->hrec[i];
            if(hr->type == BCF_HL_FMT)
            {
                for(int s = 0; s < bcf_hdr_nsamples(hdr); ++s)
                {
                    std::string header_line = "##INFO=<";
                    bool has_values = false;
                    for(int fi = 0; fi < hr->nkeys; ++fi)
                    {
                        if(has_values && fi < hr->nkeys - 1)
                        {
                            header_line += ",";
                        }
                        has_values = true;

                        if(strcmp(hr->keys[fi], "IDX") == 0)
                        {
                            continue;
                        }
                        if(strcmp(hr->keys[fi], "ID") == 0)
                        {
                            const std::string fname = hr->vals[fi];
                            const std::string sname = bcf_hdr_nsamples(hdr) > 1 ? std::string(hdr->samples[s]) : sample;
                            const std::string iname = sname + "_" + fname;
                            format_names.insert(fname);
                            header_line += "ID=";
                            header_line += iname;
                        }
                        else
                        {
                            header_line += hr->keys[fi];
                            header_line += "=";
                            header_line += hr->vals[fi];
                        }
                    }
                    header_line += ">";
                    bcf_hdr_append(output_header.get(), header_line.c_str());
                }
            }
        }
        bcf_hdr_sync(output_header.get());

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

        int nl = 1;
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

            const std::string vchr = bcfhelpers::getChrom(hdr, line);

            bcf_unpack(line, BCF_UN_ALL);

            // translate FORMAT into INFO
            for(auto const &f : format_names)
            {
                auto fmt = bcf_get_fmt(output_header.get(), line, f.c_str());

                if(fmt)
                {
                    for(int s = 0; s < bcf_hdr_nsamples(hdr); ++s)
                    {
                        const std::string sname = bcf_hdr_nsamples(hdr) > 1 ? std::string(hdr->samples[s]) : sample;
                        const std::string iname = sname + "_" + f;
                        switch(fmt->type)
                        {
                            case BCF_BT_INT8:
                            case BCF_BT_INT16:
                            case BCF_BT_INT32:
                            {
                                auto values = bcfhelpers::getFormatInts(output_header.get(),
                                                                        line, f.c_str(), s);

                                bcf_update_info_int32(output_header.get(), line, iname.c_str(), values.data(),
                                                      (int) values.size());
                                break;
                            }
                            case BCF_BT_FLOAT:
                            {
                                auto values = bcfhelpers::getFormatFloats(output_header.get(),
                                                                          line, f.c_str(), s);

                                bcf_update_info_float(output_header.get(), line, iname.c_str(), values.data(),
                                                      (int) values.size());
                                break;
                            }
                            case BCF_BT_CHAR:
                            {
                                auto value = bcfhelpers::getFormatString(output_header.get(),
                                                                         line, f.c_str(), s);

                                if(value != ".")
                                {
                                    bcf_update_info_string(output_header.get(), line, iname.c_str(), value.c_str());
                                }
                                break;
                            }
                            default:
                                error("Unsupported format field %s at %s:%i", f.c_str(), vchr.c_str(), line->pos);
                                break;
                        }
                    }
                    // if FIRST gt is kept, don't remove GT format
                    if(f != "GT" || gt_output != FIRST)
                    {
                        bcf_update_format(hdr, line, f.c_str(), NULL, 0, BCF_HT_STR);
                    }
                }
            }

            line->n_sample = 1;

            std::list<int> alleles_to_write;
            auto updateGT = [&output_header, gt_output](bcf1_t *v, bool ref)
            {
                int gtvalue = ref ? 0 : 1;
                switch(gt_output)
                {
                    case HET:
                    {
                        int gts[2] = {bcf_gt_unphased(0), bcf_gt_unphased(gtvalue)};
                        bcf_update_genotypes(output_header.get(), v, gts, 2);
                        break;
                    }
                    case HOM:
                    {
                        int gts[2] = {bcf_gt_unphased(gtvalue), bcf_gt_unphased(gtvalue)};
                        bcf_update_genotypes(output_header.get(), v, gts, 2);
                        break;
                    }
                    case HEMI:
                    {
                        int gts[1] = {bcf_gt_unphased(gtvalue)};
                        bcf_update_genotypes(output_header.get(), v, gts, 1);
                        break;
                    }
                    case HALF:
                    {
                        int gts[2] = {bcf_gt_missing, bcf_gt_unphased(gtvalue)};
                        bcf_update_genotypes(output_header.get(), v, gts, 2);
                        break;
                    }
                    default: // FIRST handled above
                        break;
                }
            };

            // dup here so we can replace the GT without segfault
            line = bcf_dup(line);
            bcf_unpack(line, BCF_UN_ALL);
            // no ALTS or FIRST genotype output (where we don't split alleles)
            if(line->n_allele == 1 || gt_output == FIRST)
            {
                updateGT(line, true);
                bcf_write(writer, output_header.get(), line);
            }
            else if(line->n_allele > 1)
            {
                updateGT(line, false);
                std::list<std::string> alleles;
                for(int a = 1; a < line->n_allele; ++a)
                {
                    alleles.emplace_back(line->d.allele[a]);
                }

                const char *alleles_to_set[2];
                const std::string ref = line->d.allele[0];
                alleles_to_set[0] = ref.c_str();
                for(auto const &al : alleles)
                {
                    alleles_to_set[1] = al.c_str();
                    bcf_update_alleles(output_header.get(), line, alleles_to_set, 2);
                    bcf_write(writer, output_header.get(), line);
                }
            }
            else
            {
                std::cerr << "Ignoring record with no alleles at " << vchr << ":" << line->pos << std::endl;
            }
            bcf_destroy(line);
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

