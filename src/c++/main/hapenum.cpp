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
 * \brief Enumerate possible Haplotypes from a VCF
 *
 * \file hapenum.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>


#include "Version.hh"
#include "Variant.hh"
#include "VariantInput.hh"
#include "Haplotype.hh"
#include "GraphReference.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"

#include <fstream>

// error needs to come after program_options.
#include "Error.hh"

using namespace variant;
using namespace haplotypes;

// specialize dot formatters
namespace graphutil
{

template <>
struct dotFormatter<ReferenceNode>
{
    void operator()(std::ostream & o, ReferenceNode const & rn)
    {
        switch(rn.type)
        {
            case ReferenceNode::source:
                o << "label=\"start: "
                  << stringutil::formatPos(rn.chr.c_str(), rn.end)
                  << "\", shape=diamond, ";
                break;
            case ReferenceNode::sink:
                o << "label=\"end: "
                  << stringutil::formatPos(rn.chr.c_str(), rn.end)
                  << "\", shape=triangle, ";
                break;
            case ReferenceNode::alternative:
                o << "label=\"alt: "
                  << rn.alt << " "
                  << stringutil::formatPos(rn.chr.c_str(), rn.start, rn.end)
                  << "\", shape=parallelogram, ";
                break;
            case ReferenceNode::homref:
                o << "label=\"homref "
                  << stringutil::formatPos(rn.chr.c_str(), rn.start, rn.end)
                  << "\", shape=box, ";
                break;
            default:
                o << "label=\"unknown "
                  << stringutil::formatPos(rn.chr.c_str(), rn.start, rn.end)
                  << "\", shape=egg, ";
                break;
        }
         o << "penwidth=2, ";
        switch(rn.color)
        {
            case ReferenceNode::red:
                o << "color=red";
                break;
            case ReferenceNode::blue:
                o << "color=blue";
                break;
            case ReferenceNode::black:
                o << "color=black";
                break;
            default:
                o << "color=cadetblue4";
                break;
        }
    }
};

} // namespace graphutil

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string ref_fasta;
    std::string file;
    std::string sample;

    std::string out_dotfile = "";
    std::string out_fasta = "";

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;

    // = max 12 unphased hets in segment
    int max_n_haplotypes = 4096;
    bool apply_filters = true;
    bool preprocess = true;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")
            ("input-file", po::value<std::string>(), "The input VCF/BCF file (use file:sample to specify a sample)")
            ("output-dot", po::value<std::string>(), "Write a dot file with the reference graph.")
            ("output-fasta", po::value<std::string>(), "Write a fasta file with all possible haplotypes.")
            ("location,l", po::value<std::string>(), "The location / subset.")
            ("reference,r", po::value<std::string>(), "The reference fasta file.")
            ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
            ("preprocess,P", po::value<bool>(), "Preprocess variants")
            ("max-n-haplotypes", po::value<int>(), "Maximum number of haplotypes to enumerate.")
        ;

        po::positional_options_description popts;
        popts.add("input-file", 1);

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
            std::cout << "hapenum version " << HAPLOTYPES_VERSION << "\n";
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
            std::vector<std::string> v;
            stringutil::split(file, v, ":");

            // in case someone passes a ":"
            assert(v.size() > 0);

            file = v[0];

            if(v.size() > 1)
            {
                sample = v[1];
            }
        }

        if (vm.count("reference"))
        {
            ref_fasta = vm["reference"].as< std::string >();
        }
        else
        {
            error("Please specify a reference file name.");
        }

        if (vm.count("output-dot"))
        {
            out_dotfile = vm["output-dot"].as< std::string >();
        }

        if (vm.count("output-fasta"))
        {
            out_fasta = vm["output-fasta"].as< std::string >();
        }

        if (vm.count("location"))
        {
            stringutil::parsePos(vm["location"].as< std::string >(), chr, start, end);
        }
        else
        {
            error("Please specify a location.");
        }

        if (vm.count("apply-filters"))
        {
            apply_filters = vm["apply-filters"].as< bool >();
        }

        if (vm.count("preprocess"))
        {
            preprocess = vm["preprocess"].as< bool >();
        }

        if (vm.count("max-n-haplotypes"))
        {
            max_n_haplotypes = vm["max-n-haplotypes"].as< int >();
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
        r.setApplyFilters(apply_filters);
        int ix = r.addSample(file.c_str(), sample.c_str());
        VariantInput vi(
            ref_fasta.c_str(),
            preprocess,           // bool leftshift
            false,           // bool refpadding
            false,          // bool trimalleles = false,
            false,          // bool splitalleles = false,
            0,              // int mergebylocation = false,
            preprocess,           // bool uniqalleles = false,
            false,          // bool calls_only = true,
            false           // bool homref_split = false
        );
        vi.getProcessor().setReader(r, VariantBufferMode::buffer_block, 100);

        std::list<Variants> vars;
        vi.get(chr.c_str(), start, end, vars);
        std::vector<ReferenceNode> nodes;
        std::vector<ReferenceEdge> edges;

        GraphReference gr(ref_fasta.c_str());
        gr.makeGraph(vars, ix, nodes, edges);

        if(out_dotfile != "")
        {
            graphutil::writeGraphDot(nodes, edges, out_dotfile.c_str());
        }

        if(out_fasta != "")
        {
            std::vector<Haplotype> haps;
            std::vector<uint64_t> nodes_used;
            std::ofstream fout(out_fasta.c_str());

            gr.enumeratePaths(chr.c_str(), start, end, 
                              nodes, edges, 
                              haps, 0, (size_t)-1, max_n_haplotypes,
                              &nodes_used);

            int i = 0;
            for(auto const & hap : haps)
            {
                std::string nodes_string = "[";

                int het_mask = 1;

                for(size_t n_id = 0; n_id < nodes.size(); ++n_id)
                {
                    // nodes_used[i] gives a bit-mask of which het nodes where used 
                    auto & n = nodes[n_id];
                    // debug-print all nodes
                    // std::cerr << n << "\n";
                    if(n.type != ReferenceNode::alternative)
                    {
                        continue;
                    }

                    bool node_was_used = false;
                    if(n.color == ReferenceNode::black)
                    {
                        node_was_used = true;
                    } 
                    else
                    {
                        node_was_used = (nodes_used[i] & het_mask) != 0;
                        het_mask <<= 1;
                    }

                    if(node_was_used)
                    {
                        if(nodes_string != "[")
                        {
                            nodes_string += " ";
                        }
                        nodes_string += std::to_string(n_id) + ":" + n.repr();                        
                    }
                }
                nodes_string += "]";

                fout << ">hap_" << i++ << ":" << stringutil::formatPos(chr, start, end)
                     << ":hb=" << stringutil::formatPos(hap.chr(), hap.start(), hap.end())
                     << " " + nodes_string
                     << "\n";
                fout << hap.seq(start, end) << "\n";
            }
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

