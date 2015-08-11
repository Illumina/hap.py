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

#include <htslib/hts.h>

#include "Error.hh"

#include "IntervalTree.h"

using namespace variant;

namespace std
{
    template<typename T, typename ...Args>
    unique_ptr<T> make_unique( Args&& ...args )
    {
        return unique_ptr<T>( new T( forward<Args>(args)... ) );
    }
}

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::vector<std::string> files;
    std::string output;
    std::string output_vcf;
    std::string ref_fasta;

    std::map<std::string, IntervalTree<std::string, int64_t> > regions;
    std::set<std::string> region_names;

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;
    int64_t rlimit = -1;

    int64_t message = -1;

    bool apply_filters = false;
    bool count_homref = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::vector< std::string> >(), "The input files")
            ("output-file,o", po::value<std::string>(), "The output file name (JSON Format).")
            ("output-vcf,v", po::value<std::string>(), "Annotated VCF file (with bed annotations).")
            ("reference,r", po::value<std::string>(), "The reference fasta file (needed only for VCF output).")
            ("location,l", po::value<std::string>(), "Start location.")
            ("regions,R", po::value< std::vector<std::string> >(),
                "Region bed file. You can attach a label by prefixing with a colon, e.g. -R FP2:false-positives-type2.bed")
            ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
            ("message-every", po::value<int64_t>(), "Print a message every N records.")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("count-homref", po::value<bool>(), "Count homref locations.")
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
            files = vm["input-file"].as< std::vector<std::string> >();
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

        if (vm.count("regions"))
        {
            std::vector<std::string> rnames = vm["regions"].as< std::vector<std::string> >();
            std::map<std::string, std::vector< Interval<std::string, int64_t > > > intervals;
            for(std::string const & f : rnames)
            {
                std::vector<std::string> v;
                stringutil::split(f, v, ":");

                std::string filename, label = "";

                // in case someone passes a ":"
                if(v.size() == 0)
                {
                    error("Invalid region name: %s", f.c_str());
                }

                if(v.size() > 1)
                {
                    label = v[0];
                    filename = v[1];
                }
                else
                {
                    filename = v[0];
                    label = boost::filesystem::path(filename).stem().string();
                }

                region_names.insert(label);

                htsFile * bedfile = NULL;

                if (stringutil::endsWith(filename, ".gz"))
                {
                    bedfile = hts_open(filename.c_str(), "rz");
                }
                else
                {
                    bedfile = hts_open(filename.c_str(), "r");
                }

                size_t icount = 0;
                kstring_t l;
                l.l = l.m = 0; l.s = NULL;
                while ( hts_getline(bedfile, 2, &l) > 0 )
                {
                    std::string line(l.s);
                    v.clear();
                    stringutil::split(line, v, "\t");
                    // we want >= 3 columns
                    if(v.size() >= 3)
                    {
                        auto chr_it = intervals.find(v[0]);
                        if (chr_it == intervals.end())
                        {
                            chr_it = intervals.emplace(
                                v[0],
                                std::vector< Interval<std::string, int64_t > >()).first;
                        }
                        // intervals are both zero-based
                        int64_t start = std::stoll(v[1]), stop = std::stoll(v[2])-1;
                        if (start > stop)
                        {
                            std::cerr << "[W] ignoring invalid interval in " << filename << " : " << line << "\n";
                        }
                        chr_it->second.push_back(Interval<std::string, int64_t>(start, stop, label));
                        ++icount;
                    }
                    else if(line != "" && line != "\n")
                    {
                        std::cerr << "[W] ignoring mis-formatted input line in " << filename << " : " << line << "\n";
                    }
                }
                free(l.s);
                hts_close(bedfile);
                std::cerr << "Added region file '" << filename << "' as '" << label << "' (" << icount << " intervals)" << "\n";
            }

            for (auto & p : intervals)
            {
                regions[p.first] = IntervalTree<std::string, int64_t> (p.second);
            }
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
        r.rewind(chr.c_str(), start);

        std::map<std::string, VariantStatistics> count_map;
        std::unique_ptr<VariantWriter> writer;

        if (output_vcf != "")
        {
            writer = std::make_unique<VariantWriter>(output_vcf.c_str(), ref_fasta.c_str());
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
                    sname = p.second + "." + std::to_string(i++);
                }
                samplenames.insert(sname);
                std::cerr << "Writing '" << p.first << ":" << p.second << "' as sample '" << sname << "'" << "\n";
                writer->addSample(sname.c_str());
            }
            writer->addHeader(r);
            writer->addHeader("##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">");
            writer->addHeader("##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">");
            writer->addHeader("##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
            writer->addHeader("##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
            writer->addHeader("##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">");
            writer->addHeader("##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">");
        }

        /** local function to count variants in all samples */
        const auto count_variants = [ref_fasta, &count_map, &samples, count_homref](std::string const & name, Variants & vars, bool update_info)
        {
            int i = 0;
            std::set<int> vtypes;
            for (auto const & s : samples)
            {
                std::string key;
                if(s.second == "")
                {
                    key = name + ":" + s.first + ":" + s.second + ":" + std::to_string(i);
                }
                else
                {
                    key = name + ":" + s.second;
                }

                auto it = count_map.find(key);
                if (it == count_map.end())
                {
                    it = count_map.emplace(key, VariantStatistics(ref_fasta.c_str(), count_homref)).first;
                }
                if(update_info)
                {
                    int * types;
                    int ntypes = 0;
                    it->second.add(vars, i, &types, &ntypes);
                    for(int j = 0; j < ntypes; ++j) {
                        vtypes.insert(types[j]);
                    }
                }
                else
                {
                    it->second.add(vars, i);
                }
                ++i;
            }
            if(update_info && !vtypes.empty())
            {
                std::string s = "VTC=";
                int j = 0;
                for(int t : vtypes)
                {
                    if(j++ > 0) { s += ","; }
                    s += VariantStatistics::type2string(t);
                }
                if (vars.info != "")
                {
                    vars.info += ";";
                }
                vars.info += s;
            }
        };

        int64_t rcount = 0;

        std::string chr = "";
        IntervalTree<std::string, int64_t> * current_chr_intervals = NULL;

        while(r.advance())
        {
            Variants & v = r.current();

            if(rlimit != -1)
            {
                if(rcount >= rlimit)
                {
                    break;
                }
            }

            if(end != -1 && v.pos > end)
            {
                break;
            }

            if (!regions.empty() && (current_chr_intervals == NULL || v.chr != chr))
            {
                auto qit = regions.find(v.chr);
                if (qit != regions.end())
                {
                    current_chr_intervals = &qit->second;
                }
            }
            chr = v.chr;
            count_variants("all", v, true);
            // resolve tags
            std::set<std::string> to_count;
            std::string tag_string = "";
            if (current_chr_intervals != NULL)
            {
                std::vector< Interval<std::string, int64_t > > overlapping;
                current_chr_intervals->findOverlapping(v.pos, v.pos + v.len - 1, overlapping);
                for (auto const & iv : overlapping)
                {
                    to_count.insert(iv.value);
                }
                for (std::string const & r : to_count)
                {
                    if (tag_string != "")
                    {
                        tag_string += ",";
                    }
                    tag_string += r;
                }
            }

            std::string type = ".";
            std::string kind = ".";
            std::string gtt1 = ".";
            std::string gtt2 = ".";
            bool hapmatch = false;

            std::vector<std::string> infs;
            stringutil::split(v.info, infs, ";");
            v.info = "";
            for (std::string & i : infs)
            {
                if(infs.size() < 2)
                {
                    continue;
                }
                std::vector<std::string> ifo;
                stringutil::split(i, ifo, "=");
                if(ifo[0] == "HapMatch") {
                    hapmatch = true;
                } else if(ifo.size() < 2) {
                    v.info += i + ";";
                } else if(ifo[0] == "type") {
                    type = ifo[1];
                } else if(ifo[0] == "kind") {
                    kind = ifo[1];
                } else if(ifo[0] == "gtt1") {
                    gtt1 = ifo[1];
                } else if(ifo[0] == "gtt2") {
                    gtt2 = ifo[1];
                } else if(ifo[0] != "ctype") {
                    v.info += i + ";";
                }
            }

            if(hapmatch && type != "TP") {
                kind = "hapmatch__" + type + "__" + kind;
                type = "TP";
            }

            if(!regions.empty() && tag_string.empty() && type == "FP" && kind == "missing") {
                type = "UNK";
            }

            v.info += std::string("type=") + type;
            v.info += std::string(";kind=") + kind;
            v.info += std::string(";gtt1=") + gtt1;
            v.info += std::string(";gtt2=") + gtt2;
            if(!tag_string.empty()) {
                v.info += std::string(";Regions=") + tag_string;
            }

            if (writer)
            {
                writer->put(v);
            }

            count_variants(type + ":" + kind + ":" + tag_string, v, false);

            if(message > 0 && (rcount % message) == 0)
            {
                std::cout << stringutil::formatPos(v.chr.c_str(), v.pos) << ": " << v << "\n";
            }
            ++rcount;
        }

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
