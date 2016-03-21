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
 *  \brief Read a spacer-separated file and compute a ROC curve
 *
 *
 * \file roc.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>

#include "Version.hh"
#include "helpers/StringUtil.hh"

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>

#include "Error.hh"

enum class tag_t { tp, fp, fn };

struct RocData {
    void add(double value, tag_t tag)
    {
        if(reverse) value = -value;
        data.push_back(std::make_pair(value, tag));
    }

    void print_roc(std::ostream & out) {
        if(data.empty())
        {
            return;
        }
        // pair gets sorted by first element, then second
        std::sort(data.begin(), data.end());
        int total_tp = 0;
        int total_fp = 0;
        int total_fn = 0;
        for (auto r : data) {
            switch(r.second)
            {
                case tag_t::tp: ++total_tp; break;
                case tag_t::fp: ++total_fp; break;
                case tag_t::fn: ++total_fn; break;
                default: break;
            }
        }

        out << label << "\t" << "tp" << "\t" << "fp" << "\t" << "fn" << "\t" << "precision" << "\t" << "recall" << "\n";
        double current_value = data[0].first - 1;
        double tp = total_tp;
        double fp = total_fp;
        double fn = total_fn;
        for (auto r : data) {
            if(r.first > current_value)
            {
                current_value = r.first;
                double precision = 1.0;
                double recall = 0;
                if(tp + fp > 0)
                {
                    precision = tp / (tp + fp);
                }
                if(tp + fn > 0)
                {
                    recall = tp / (tp + fn);
                }

                out << (reverse ? -current_value : current_value) << "\t" << tp << "\t" << fp << "\t" << fn << "\t" << precision << "\t" << recall << "\n";
            }

            switch(r.second)
            {
                case tag_t::fn:
                    break;
                case tag_t::tp:
                    --tp;
                    ++fn;
                    break;
                case tag_t::fp:
                    --fp;
                    break;
            }
        }
    }

    std::vector<std::pair<double, tag_t> > data;
    std::string label = "value";

    bool reverse = false;
};


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::vector<std::string> files;
    std::string output = "-";
    std::string sep = "\t";
    int header_lines = 1;
    std::string value = "";
    int value_column = 0;
    std::string tag = "";
    int tag_column = 1;
    std::string filter = "";
    int filter_column = -1;
    std::string filter_name = "";

    bool verbose = false;
    bool reverse = false;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "show version")
            ("verbose", "show verbose information (to stderr)")
            ("input-file", po::value<std::vector< std::string> >(), "The input files")
            ("output-file,o", po::value<std::string>(), "Output file name, defaults to - / write to stdout")
            ("separator,s", po::value<std::string>(), "separator character (default: '\\t' for reading tsv)")
            ("header-lines,H", po::value<int>(), "lines to skip before starting to read")
            ("value,v", po::value<std::string>(), "value column name")
            ("value-column", po::value<int>(), "value column number")
            ("reverse,R", po::value<bool>(), "Reverse counting for score (default: higher scores are better)")
            ("tag,t", po::value<std::string>(), "tag column name")
            ("tag-column", po::value<int>(), "tag column number. Tags must be TP/FP/FN, lines with different tags will be ignored")
            ("filter,f", po::value<std::string>(), "filter column name")
            ("filter-column", po::value<int>(), "filter column number. This is used if we the value we are varying is a threshold for a certain filter.")
            ("filter-name,n", po::value<std::string>(), "filter name if value is threshold for this filter")
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
            std::cout << "roc version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("verbose"))
        {
            std::cerr << "Verbose mode enabled" << "\n";
            verbose = true;
        }

        if (vm.count("input-file"))
        {
            files = vm["input-file"].as< std::vector<std::string> >();
        }

        if (vm.count("output-file"))
        {
            output = vm["output-file"].as<std::string>();
        }

        if (vm.count("separator"))
        {
            sep = vm["separator"].as< std::string >();
        }

        if (vm.count("header-lines"))
        {
            header_lines = vm["header-lines"].as<int>();
        }

        if (vm.count("value-column"))
        {
            value_column = vm["value-column"].as<int>() - 1;
        }
        if (vm.count("value"))
        {
            value = vm["value"].as<std::string>();
        }
        if (vm.count("reverse"))
        {
            reverse = vm["reverse"].as<bool>();
        }

        if (vm.count("tag-column"))
        {
            tag_column = vm["tag-column"].as<int>() - 1;
        }
        if (vm.count("tag"))
        {
            tag = vm["tag"].as<std::string>();
        }

        if (vm.count("filter-column"))
        {
            filter_column = vm["filter-column"].as<int>() - 1;
        }
        if (vm.count("filter"))
        {
            filter = vm["filter"].as<std::string>();
        }
        if (vm.count("filter-name"))
        {
            filter_name = vm["filter-name"].as<std::string>();
        }
    }
    catch (po::error & e)
    {
        std::cerr << e.what() << "\n";
        return 1;
    }


    RocData data;
    data.reverse = reverse;
    if(!value.empty())
    {
        data.label = value;
    }

    if(files.size() == 0) {
        files.push_back("-");
    }

    int total_tp = 0;
    int total_fp = 0;
    int total_fn = 0;
    int total_filtered = 0;
    int total_ignored = 0;
    std::map<std::string, int> filter_stats;
    std::map<std::string, int> filter_stats_tp;
    std::map<std::string, int> filter_stats_fp;
    std::map<std::string, int> filter_stats_fn;

    for (auto i : files) {
        std::istream * in = NULL;
        bool deleteme = false;
        if(i == "-") {
            in = &std::cin;
        } else {
            in = new std::ifstream(i.c_str());
            deleteme = true;
        }

        // read input data
        int min_columns = std::max(value_column, tag_column);
        min_columns = std::max(filter_column, min_columns);

        int hc = header_lines;
        while(in->good()) {
            std::string line;
            std::getline(*in, line);
            std::vector<std::string> v;
            stringutil::split(line, v, sep, true);

            if(hc > 0) {
                if (!value.empty() || !tag.empty() || !filter.empty())
                {
                    for (size_t i = 0; i < v.size(); ++i)
                    {
                        v[i] = stringutil::replaceAll(v[i], "\r", "");
                        if(!v[i].empty() && v[i] == value) {
                            value_column = i;
                        }
                        if(!v[i].empty() && v[i] == tag) {
                            tag_column = i;
                        }
                        if(!v[i].empty() && v[i] == filter) {
                            filter_column = i;
                        }
                    }
                }
                min_columns = std::max(value_column, tag_column);
                min_columns = std::max(filter_column, min_columns);
                --hc;
                continue;
            }

            if(((int)v.size()) < min_columns + 1) {
                ++total_ignored;
                continue;
            }

            std::string & ltag = v[tag_column];
            boost::algorithm::to_lower(ltag);

            // true if filtered by other filter than the one we're looking at. These go to the beginning.
            bool filtered_other = false;

            if (filter_column >= 0)
            {
                std::vector<std::string> filters;
                std::string fcol = v[filter_column];
                if(fcol == "." || fcol == "PASS") {
                    fcol = "";
                }
                stringutil::split(fcol, filters, ";,");
                for (auto f : filters)
                {
                    auto q = filter_stats.find(f);
                    if(q != filter_stats.end()) {
                        q->second++;
                    } else {
                        filter_stats[f] = 1;
                    }
                    if (ltag == "tp")
                    {
                        auto q = filter_stats_tp.find(f);
                        if(q != filter_stats_tp.end()) {
                            q->second++;
                        } else {
                            filter_stats_tp[f] = 1;
                        }
                    } else if (ltag == "fp") {
                        auto q = filter_stats_fp.find(f);
                        if(q != filter_stats_fp.end()) {
                            q->second++;
                        } else {
                            filter_stats_fp[f] = 1;
                        }
                    } else if (ltag == "fn") {
                        auto q = filter_stats_fn.find(f);
                        if(q != filter_stats_fn.end()) {
                            q->second++;
                        } else {
                            filter_stats_fn[f] = 1;
                        }
                    }
                }

                if(!filter_name.empty())
                {
                    auto it = filters.begin();
                    std::vector<std::string> _filters_to_remove;
                    stringutil::split(filter_name, _filters_to_remove, ";,");
                    std::set<std::string> filters_to_remove;
                    for (auto f : _filters_to_remove) {
                        filters_to_remove.insert(f);
                        if(f == "*")
                        {
                            filters.clear();
                            break;
                        }
                    }
                    while(it != filters.end())
                    {
                        if(filters_to_remove.count(*it) > 0) {
                            filters.erase(it);
                            it = filters.begin();
                        } else {
                            ++it;
                        }
                    }
                }
                if (!filters.empty())
                {
                    ++total_filtered;
                    filtered_other = true;
                }
            }

            tag_t ttag = tag_t::fn;

            if(ltag.substr(0, 2) == "tp") {
                ttag = tag_t::tp;
                ++total_tp;
            } else if(ltag.substr(0, 2) == "fp") {
                ttag = tag_t::fp;
                ++total_fp;
            } else if(ltag.substr(0, 2) == "fn") {
                ttag = tag_t::fn;
                ++total_fn;
            } else {
                ++total_ignored;
                continue;
            }

            double xvalue = atof(v[value_column].c_str());
            if(v[value_column] == "")
            {
                xvalue = 0;
            }
            if (filtered_other)
            {
                xvalue = std::numeric_limits<double>::min();
            }
            data.add(xvalue, ttag);
        }
        if (deleteme)
        {
            delete in;
        }
    }
    if(verbose) {
        std::cerr << "tp: " << total_tp << " fp: " << total_fp << " fn: " << total_fn
                  << " filtered: " << total_filtered << " ignored: " << total_ignored << "\n";
        for (auto i : filter_stats) {
            int f_tp = 0;
            if (filter_stats_tp.find(i.first) != filter_stats_tp.end())
            {
                f_tp = filter_stats_tp[i.first];
            }
            int f_fp = 0;
            if (filter_stats_fp.find(i.first) != filter_stats_fp.end())
            {
                f_fp = filter_stats_fp[i.first];
            }
            std::cerr << "F:" << i.first << " : " << i.second
                      << " -- tp: " << f_tp
                      << " fp: " << f_fp
                      << "\n";
        }
    }
    if(output == "-") {
        data.print_roc(std::cout);
    } else {
        std::ofstream o(output);
        data.print_roc(o);
    }

    return 0;
}
