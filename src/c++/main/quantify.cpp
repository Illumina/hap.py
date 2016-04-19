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
#include "QuantifyRegions.hh"

using namespace variant;


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace bf = boost::filesystem;

    std::string file;
    std::string output_roc;
    std::string output_vcf;
    std::string ref;
    std::string only_regions;
    std::string qtype = "xcmp";
    std::string qq = "QUAL";
    std::string roc_filter = "";
    double roc_delta = 0.1;

    // limits
    std::string chr;
    int64_t start = -1;
    int64_t end = -1;
    int64_t rlimit = -1;

    int64_t message = -1;

    bool apply_filters = false;
    bool count_homref = false;
    bool output_vtc = true;
    bool clean_info = true;
    bool fixchr = false;
    bool output_rocs = true;

    int threads = 1;
    int blocksize = 20000;

    try
    {
        QuantifyRegions regions;

        try
        {
            // Declare the supported options.
            po::options_description desc("Allowed options");
            desc.add_options()
                ("help,h", "produce help message")
                ("version", "Show version")
                ("input-file", po::value<std::string>(), "The input file")
                ("output-file,o", po::value<std::string>(), "The output file name (TSV Format).")
                ("output-vcf,v", po::value<std::string>(), "Annotated VCF file (with bed annotations).")
                ("output-summary", po::value<std::string>(), "Output a summary table with TP / FP / FN / UNK counts, precision, recall, etc.")
                ("output-filter-rocs", po::value<bool>(), "Output ROC levels for filters.")
                ("roc-filter", po::value<std::string>(), "Ignore certain filters when creating a ROC.")
                ("roc-delta", po::value<double>(), "Minium spacing of levels on ROC QQ trace.")
                ("qq", po::value<std::string>(), "Field to use for QQ (ROC quantity). Can be QUAL / GQ / ... / any INFO field name.")
                ("reference,r", po::value<std::string>(), "The reference fasta file (needed only for VCF output).")
                ("location,l", po::value<std::string>(), "Start location.")
                ("regions,R", po::value< std::vector<std::string> >(),
                    "Region bed file. You can attach a label by prefixing with a colon, e.g. -R FP2:false-positives-type2.bed")
                ("type", po::value<std::string>(), "Quantification method to use. Current choices are xcmp or ga4gh.")
                ("only,O", po::value< std::string >(), "Bed file of locations (equivalent to -R in bcftools)")
                ("limit-records", po::value<int64_t>(), "Maximum umber of records to process")
                ("message-every", po::value<int64_t>(), "Print a message every N records.")
                ("apply-filters,f", po::value<bool>(), "Apply filtering in VCF.")
                ("count-homref", po::value<bool>(), "Count homref locations.")
                ("output-vtc", po::value<bool>(), "Output variant types counted (debugging).")
                ("clean-info", po::value<bool>(), "Set to zero to preserve INFO fields (default is 1)")
                ("output-rocs", po::value<bool>(), "Output ROCs with full set of levels ofver QQ values (default is 1, disable for more concise output)")
                ("fix-chr-regions", po::value<bool>(), "Add chr prefix to regions if necessary (default is off).")
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
                output_roc = vm["output-file"].as< std::string >();
            }

            if (vm.count("output-vcf"))
            {
                output_vcf = vm["output-vcf"].as< std::string >();
            }

            if (vm.count("reference"))
            {
                ref = vm["reference"].as< std::string >();
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

            if (vm.count("type"))
            {
                qtype = vm["type"].as< std::string >();
            }

            if (vm.count("qq"))
            {
                qq = vm["qq"].as< std::string >();
            }

            if (vm.count("roc-filter"))
            {
                roc_filter = vm["roc-filter"].as< std::string >();
            }

            if (vm.count("roc-delta"))
            {
                roc_delta = vm["roc-delta"].as< double >();
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

            if (vm.count("clean-info"))
            {
                clean_info = vm["clean-info"].as< bool >();
            }

            if (vm.count("output-rocs"))
            {
                output_rocs = vm["output-rocs"].as< bool >();
            }

            if (vm.count("fix-chr-regions"))
            {
                fixchr = vm["fix-chr-regions"].as< bool >();
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

            if (output_roc == "")
            {
                std::cerr << "Please specify an output file.\n";
                return 1;
            }

            if (vm.count("regions"))
            {
                std::vector<std::string> rnames = vm["regions"].as< std::vector<std::string> >();
                regions.load(rnames, fixchr);
            }

        }
        catch (po::error & e)
        {
            std::cerr << e.what() << "\n";
            return 1;
        }

        FastaFile ref_fasta(ref.c_str());
        bcf_srs_t * reader = bcf_sr_init();
        reader->require_index = 1;
        reader->collapse = COLLAPSE_NONE;
        reader->streaming = 0;
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
                std::cerr << "starting at " << chr << ":" << start << "\n";
            }
            if(success < 0)
            {
                error("Cannot seek to %s:%i", chr.c_str(), start);
            }
        }

        bcf_hdr_t * hdr = reader->readers[0].header;

        std::string qparams = "";
        if(output_vtc)
        {
            qparams += "output_vtc;";
        }
        if(regions.hasRegions())
        {
            qparams += "count_unk;";
        }
        if(count_homref)
        {
            qparams += "count_homref;";
        }
        if(clean_info)
        {
            qparams += "clean_info;";
        }
        if(!qq.empty() || qq == ".")
        {
            qparams += "QQ:" + qq + ";";
        }
        qparams += "extended_counts;";

        std::unique_ptr<BlockQuantify> p_bq(std::move(makeQuantifier(hdr, ref_fasta, qtype, qparams)));
        p_bq->rocFiltering(roc_filter);

        // update the header
        p_bq->updateHeader(hdr);

        htsFile * writer = nullptr;

        if (output_vcf != "")
        {

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
        /** async stuff. each block can be counted in parallel, but we need to
         *  write out the variants sequentially.
         *  Therefore, we keep a future for each block to be able to join
         *  when it's processed
         */
        std::queue<std::pair <
            std::future<void>,
            std::unique_ptr<BlockQuantify>
        >> blocks;


        int truth_sample_id = -1;
        int query_sample_id = -1;

        if(bcf_hdr_nsamples(hdr) < 2)
        {
            if(!output_roc.empty())
            {
                std::cerr << "[W] not enough samples in input file. Switching off ROC table calculation" << "\n";
            }
            output_roc = "";
        }
        else
        {
            for(int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
            {
                if(strcmp(hdr->samples[i], "TRUTH") == 0)
                {
                    truth_sample_id = i;
                }
                if(strcmp(hdr->samples[i], "QUERY") == 0)
                {
                    query_sample_id = i;
                }
            }
            if(truth_sample_id < 0 || query_sample_id < 0)
            {
                if(!output_roc.empty())
                {
                    std::cerr << "[W] Input VCF does not have TRUTH and QUERY samples. Switching off ROC table calculation" << "\n";
                }
                output_roc = "";
            }
        }

        std::unique_ptr<std::map<std::string, roc::Roc> > roc_map;
        if(!output_roc.empty())
        {
            roc_map.reset(new std::map<std::string, roc::Roc>());
        }

        /** this is where things actually get written to files */
        auto output_counts = [&writer, &roc_map, &blocks, hdr, truth_sample_id, query_sample_id](int min_size) {
            while(blocks.size() > (unsigned )min_size)
            {
                // make sure we have run this block
                blocks.front().first.get();

                // output variants
                if(writer)
                {
                    auto const & variants = blocks.front().second->getVariants();
                    for(auto & v : variants)
                    {
                        bcf_write1(writer, hdr, v);
                    }
                }

                // update ROC data
                if(roc_map)
                {
                    auto const & rm = blocks.front().second->getRocs();
                    for(auto const & r : rm)
                    {
                        auto it = roc_map->find(r.first);
                        if (it == roc_map->end()) {
                            roc_map->insert(r);
                        }
                        else
                        {
                            it->second.add(r.second);
                        }
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

            bcf_unpack(line, BCF_UN_INFO);
            regions.annotate(hdr, line);

            current_chr = vchr;

            p_bq->add(bcf_dup(line));

            ++vars_in_block;
            if(vars_in_block > blocksize)
            {
                std::future<void> f = std::async(std::launch::async, &BlockQuantify::count, p_bq.get());
                // clear / write out some blocks (make sure we have at least 2xthreads tasks left)
                output_counts(threads);
                blocks.emplace(std::move(f), std::move(p_bq));
                p_bq = std::move(makeQuantifier(hdr, ref_fasta, qtype, qparams));
                p_bq->rocFiltering(roc_filter);
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

        if(writer)
        {
            hts_close(writer);
        }
        bcf_sr_destroy(reader);

        if(roc_map)
        {
            auto makeHeader = [qq](std::ostream & out)
            {
                out << "Type";
                out << "\t" << "Subtype";
                out << "\t" << "Genotype";
                out << "\t" << "Filter";
                out << "\t" << "Subset";
                out << "\t" << qq;
                for(int j = 0; j < roc::NDecisionTypes; ++j)
                {
                    if(j == roc::to_underlying(roc::DecisionType::FN2))
                    {
                        continue;
                    }
                    out << "\t" << roc::DecisionTypes[j];
                }

                out << "\t" << "FP.gt";
                out << "\t" << "FP.al";

                out << "\t" << "METRIC.Recall";
                out << "\t" << "METRIC.Precision";
                out << "\t" << "METRIC.F1_Score";
                out << "\t" << "METRIC.Frac_NA";
                out << "\t" << "TRUTH.TOTAL";
                out << "\t" << "QUERY.TOTAL";
                out << "\n";
            };

            std::ofstream out_roc(output_roc);
            makeHeader(out_roc);
            for(auto & m : *roc_map)
            {
                // filter ROCs just give numbers of TPs / FPs / UNKs filtered by
                // each filter at each level

                std::vector<std::string> subs;
                stringutil::split(m.first, subs, ":");
                if(subs.empty())
                {
                    // only happens when an empty string is passed
                    continue;
                }
                std::string type = "unknown";
                std::string subset = "unknown";
                std::string filter = "unknown";

                bool counts_only = false;
                if (subs.size() != 3)
                {
                    type = m.first;
                }
                else
                {
                    type = subs[1];
                    filter = subs[2];
                    if(subs[0] == "f")
                    {
                        counts_only = true;
                        subset = "*";
                    }
                    else if(subs[0].substr(0, 2) == "s|")
                    {
                        subset = subs[0].substr(2);
                    }
                    else if(subs[0] == "a")
                    {
                        subset = "*";
                    }
                    else
                    {
                        // unknown count, should not happen
                        continue;
                    }
                }

                std::list< std::pair<std::string, uint64_t> > subtypes;

                std::list <std::pair<std::string, uint64_t> > gts = {
                    {"*", 0},
                    {"het", roc::OBS_FLAG_HET},
                    {"hetalt", roc::OBS_FLAG_HETALT},
                    {"homalt", roc::OBS_FLAG_HOMALT},
                };

                std::list <std::pair<std::string, uint64_t> > sts = {
                    {"*", 0},
                };
                // all types, all genotypes

                // SNP: TI / TV
                if(type == "SNP")
                {
                    sts.emplace_back("ti", roc::OBS_FLAG_TI);
                    sts.emplace_back("tv", roc::OBS_FLAG_TV);
                }
                else if(type == "INDEL")
                {
                    sts.emplace_back("I1_5", roc::OBS_FLAG_I1_5);
                    sts.emplace_back("I6_15", roc::OBS_FLAG_I6_15);
                    sts.emplace_back("I16_PLUS", roc::OBS_FLAG_I16_PLUS);
                    sts.emplace_back("D1_5", roc::OBS_FLAG_D1_5);
                    sts.emplace_back("D6_15", roc::OBS_FLAG_D6_15);
                    sts.emplace_back("D16_PLUS", roc::OBS_FLAG_D16_PLUS);
                    sts.emplace_back("C1_5", roc::OBS_FLAG_C1_5);
                    sts.emplace_back("C6_15", roc::OBS_FLAG_C6_15);
                    sts.emplace_back("C16_PLUS", roc::OBS_FLAG_C16_PLUS);
                }

                for(auto const & st : sts)
                {
                    for(auto const & gt : gts)
                    {
                        subtypes.emplace_back(st.first + "\t" + gt.first, st.second | gt.second);
                    }
                }

                for(auto const & st : subtypes)
                {
                    const roc::Level l = m.second.getTotals(st.second);
                    out_roc << type << "\t" << st.first << "\t" << filter << "\t" << subset << "\t";
                    l.dumpTSV(out_roc, counts_only);
                    out_roc << "\n";

                    // don't write ROCs for filters
                    // we could make this optional, but not sure it is really necessary
                    if((!counts_only) && output_rocs)
                    {
                        std::vector<roc::Level> levels;
                        m.second.getLevels(levels, roc_delta, st.second);
                        for (auto const &l2 : levels)
                        {
                            out_roc << type << "\t" << st.first << "\t" << filter << "\t" << subset << "\t";
                            l2.dumpTSV(out_roc, counts_only);
                            out_roc << "\n";
                        }
                    }
                }
            }
            out_roc.close();
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
