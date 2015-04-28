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
 * \brief use HTSLib to turn a VCF header into JSON
 *
 * \file vcfhdr2json.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

#include "json/json.h"

#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "Version.hh"
#include "Error.hh"

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string file;
    std::string output;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version", "Show version")            
            ("input-file", po::value< std::string >(), "The input files")
            ("output-file", po::value<std::string>(), "The output file name.")
        ;

        po::positional_options_description popts;
        popts.add("input-file", 1);
        popts.add("output-file", 1);

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
            std::cout << "multimerge version " << HAPLOTYPES_VERSION << "\n";
            return 0;
        }

        if (vm.count("help")) 
        {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("input-file"))
        {
            file = vm["input-file"].as< std::string > ();
        }

        if (vm.count("output-file"))
        {
            output = vm["output-file"].as< std::string >();
        }

        if(file.size() == 0)
        {
            std::cerr << "Please specify an input file.\n";
            return 1;
        }

        if (output == "")
        {
            std::cerr << "Please specify an output file.\n";
            return 1; 
        }
    } 
    catch (po::error & e)
    {
        std::cerr << e.what() << "\n";
        return 1;
    }

    try
    {
        Json::StyledWriter writer;
        htsFile * fp = bcf_open(file.c_str(), "r");
        bcf_hdr_t * hdr = bcf_hdr_read(fp);

        Json::Value root;
        Json::Value a;
        for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
        {
            a.append(hdr->samples[i]);
        }
        root["samples"] = a;

        Json::Value fields;
        for (int i = 0; i < hdr->nhrec; i++)
        {
            Json::Value field;
            field["key"] = hdr->hrec[i]->key;
            if (!hdr->hrec[i]->value)
            {
                Json::Value values;

                for (int j = 0; j < hdr->hrec[i]->nkeys; j++)
                {
                    values[hdr->hrec[i]->keys[j]] = hdr->hrec[i]->vals[j];
                }
                field["values"] = values;
            }
            else
            {
                field["value"] = hdr->hrec[i]->value;
            }
            fields.append(field);
        }
        root["fields"] = fields;

        tbx_t * tbx_idx = tbx_index_load(file.c_str());
        if ( !tbx_idx )
        {
            root["tabix"] = Json::Value::null;
        }
        else
        {
            root["tabix"] = Json::Value();
            root["tabix"]["chromosomes"] = Json::Value();

            int count = 0;
            const char ** tbx_names = tbx_seqnames(tbx_idx, &count);

            for (int i = 0; i < count; ++i)
            {
                root["tabix"]["chromosomes"].append(tbx_names[i]);
            }

            tbx_destroy(tbx_idx);
        }


        std::ofstream out(output.c_str());
        out << writer.write(root);

        bcf_close(fp);
        bcf_hdr_destroy(hdr);
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
