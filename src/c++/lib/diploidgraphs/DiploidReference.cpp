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
 * Partially phased diploid reference representation
 *
 * \file DiploidReference.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "DiploidReference.hh"

#include "GraphReference.hh"

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <set>
#include <map>

#include "Error.hh"

// #define _DEBUG_DIPLOIDREFERENCE

namespace haplotypes
{

std::ostream & operator<<(std::ostream & o, DiploidRef const & r)
{
    if(r.het)
    {
        if(r.homref)
        {
            o << "het(" << r.h1 << "|" << r.h2 << ")";
        }
        else
        {
            o << "hetalt(" << r.h1 << "|" << r.h2 << ")";
        }
    }
    else
    {
        if(r.homref)
        {
            o << "homref(" << r.h1 << ")";
        }
        else
        {
            o << "hom(" << r.h1 << ")";
        }
    }
    return o;
}

struct DiploidReferenceImpl
{
    DiploidReferenceImpl(const char * _vcfname, 
                         const char * _sample,
                         const char * _ref_fasta) : 
        vcfname(_vcfname ? _vcfname : ""), sample(_sample ? _sample : ""), ref_fasta(_ref_fasta),
        chr(""), start(-1), end(-1), max_n_paths(-1),
        gr(_vcfname, _sample, _ref_fasta)
        {}
    
    /* copyables */
    std::string vcfname;
    std::string sample;
    std::string ref_fasta;

    std::string chr;
    int64_t start;
    int64_t end;

    int max_n_paths;

    /* created by DiploidReference */
    GraphReference gr;

    // created by setRegion
    std::list< DiploidRef > di_haps;
    std::list< DiploidRef > :: iterator pos;
};

DiploidReference::DiploidReference(GraphReference const & gr) : 
    _impl(new DiploidReferenceImpl(gr.getVCFName().c_str(), gr.getSample().c_str(), gr.getFastaName().c_str()))
{
    _impl->gr.setApplyFilters(gr.getApplyFilters());
}


DiploidReference::DiploidReference(
    const char * vcfname, 
    const char * sample,
    const char * ref_fasta) : _impl(new DiploidReferenceImpl(vcfname, sample, ref_fasta))
{
}

DiploidReference::DiploidReference(DiploidReference const & rhs)
{

    _impl = new DiploidReferenceImpl(rhs._impl->vcfname.c_str(), 
                                     rhs._impl->sample.c_str(), 
                                     rhs._impl->ref_fasta.c_str());  

    _impl->gr.setApplyFilters(rhs._impl->gr.getApplyFilters());
}

DiploidReference const & DiploidReference::operator=(DiploidReference const & rhs)
{
    if(&rhs == this)
    {
        return *this;
    }
    delete _impl;

    DiploidReferenceImpl * t_impl = new DiploidReferenceImpl(rhs._impl->vcfname.c_str(), 
                                                             rhs._impl->sample.c_str(), 
                                                             rhs._impl->ref_fasta.c_str());

    _impl = t_impl;

    _impl->gr.setApplyFilters(rhs._impl->gr.getApplyFilters());
    return *this;
}

DiploidReference::~DiploidReference()
{
    delete _impl;
}

/**
 * Get reference fasta file
 */
FastaFile & DiploidReference::getRefFasta()
{
    return _impl->gr.getRefFasta();
}

FastaFile const & DiploidReference::getRefFasta() const
{
    return _impl->gr.getRefFasta();
}

/**
 * Set the maximum number of paths to enumerate from the Graph reference
 */
void DiploidReference::setNPaths(int max_n_paths)
{
    _impl->max_n_paths = max_n_paths;
}

/**
 * Set region and restart enumeration
 */
void DiploidReference::setRegion(
    const char * chr,
    int64_t start,
    int64_t end,
    std::list<variant::Variants> const * vars,
    int sample_ix)
{
    _impl->chr = chr;
    _impl->start = start;
    _impl->end = end;
    _impl->di_haps.clear();

    if(std::string(chr) != "" && start >= 0 && end >= 0 && end - start + 1 > 0)
    {
        std::vector<ReferenceNode> nodes;
        std::vector<ReferenceEdge> edges;

        size_t nhets = 0;
        if (vars)
        {
            _impl->gr.makeGraph(*vars, sample_ix, nodes, edges, &nhets);
        }
        else
        {
            _impl->gr.makeGraph(chr, start, end, nodes, edges, &nhets);
        }

#ifdef _DEBUG_DIPLOIDREFERENCE
        std::cerr << "nhets : " << nhets << " mp: " << _impl->max_n_paths << "\n";
#endif

        if(_impl->max_n_paths > 0 && pow(2.0, (double)nhets) > _impl->max_n_paths)
        {
            error("Too many het nodes (%i) at %s:%i-%i", nhets, chr, start, end);
        } 

        std::string refsq = getRefFasta().query(chr, start, end);

        std::vector<Haplotype> target;
        std::vector<std::string> nodes_used;

        _impl->gr.enumeratePaths(chr, start, end, 
                          nodes, edges, target, 
                          0, _impl->max_n_paths, -1, 
                          &nodes_used);

#ifdef _DEBUG_DIPLOIDREFERENCE
        std::cerr << "Nodes: " << "\n";
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            std::cerr << i << ": " << nodes[i] << " het:" << nodes[i].het << " filtered: " << nodes[i].filtered << "\n";
        }
        std::cerr << "Haplotypes seen: " << "\n";
        for (size_t i = 0; i < target.size(); ++i)
        {
            std::cerr << target[i].repr(start, end) << "\n";
            std::cerr << "NU: " << nodes_used[i] << "\n";
        }
#endif
        std::list<size_t> het_nodes;
        std::list<size_t> hom_nodes;
        for(size_t n = 0; n < nodes.size(); ++n)
        {
            if (nodes[n].type == ReferenceNode::alternative) 
            {
                if(   nodes[n].het
                   || nodes[n].color == ReferenceNode::red 
                   || nodes[n].color == ReferenceNode::blue )
                {
                    het_nodes.push_back(n);
                }
                else if(!nodes[n].het)
                {
                    hom_nodes.push_back(n);
                }
            }
        }

        // if we have het nodes, each pair of haplotypes must cover them all
        if(het_nodes.size() != 0)
        {
            if(_impl->max_n_paths > 0 && pow(2, het_nodes.size()) > _impl->max_n_paths)
            {
                error("Too many het nodes (%i) at %s:%i-%i", het_nodes.size(), chr, start, end);
            } 
            for (size_t p1 = 0; p1 < nodes_used.size(); ++p1)
            {
                for (size_t p2 = p1; p2 < nodes_used.size(); ++p2) // we allow matching up p1 with p1 to allow creation of hom blocks via filtered variants
                {
                    size_t hn_found = 0;
                    size_t fn_used = 0;
                    size_t fn_unused = 0;
                    for(size_t n : het_nodes)
                    {
                        if(nodes[n].filtered)
                        {
                            if(nodes_used[p1][nodes.size() - 1 - n] == '1' || nodes_used[p2][nodes.size() - 1 - n] == '1')
                            {
                                ++fn_used;
                            }
                            if(nodes_used[p1][nodes.size() - 1 - n] == '0' && nodes_used[p2][nodes.size() - 1 - n] == '0')
                            {
                                ++fn_unused;
                            }
                        }

                        if(nodes_used[p1][nodes.size() - 1 - n] + nodes_used[p2][nodes.size() - 1 - n] == '0' + '1')
                        {
                            ++hn_found;
                        }
#ifdef _DEBUG_DIPLOIDREFERENCE
                        std::cerr << "n: " << n << " p1:" << p1 << "(" << nodes_used[p1][nodes.size() - 1 - n] << ")"
                                  << " p2:" << p2 << "(" <<  nodes_used[p2][nodes.size() - 1 - n] << ")" 
                                  << " filtered: " << nodes[n].filtered << " nf: " << fn_used << " nfu: " << fn_unused
                                  << "\n";
#endif
                    }

                    // make sure hom nodes are covered as hom
                    bool hom_conflict = false;
                    for(size_t n : hom_nodes)
                    {
                        if (nodes[n].filtered)
                        {
                            // for filtered variants, we need them to either be present or gone!
                            if(nodes_used[p1][nodes.size() - 1 - n] != nodes_used[p2][nodes.size() - 1 - n])
                            {
                                hom_conflict = true;
                            }
                        }
                        else
                        {
                            // unfiltered hom variants must be present in every path
                            if(nodes_used[p1][nodes.size() - 1 - n] != nodes_used[p2][nodes.size() - 1 - n] ||
                               nodes_used[p2][nodes.size() - 1 - n] != '1')
                            {
                                hom_conflict = true;
                            }
                        }
                    }

#ifdef _DEBUG_DIPLOIDREFERENCE
                    std::cerr << "p1: " << p1 << "/" << nodes_used[p1] << " "
                              << "p2: " << p2 << "/" << nodes_used[p2] << " "
                              << "\n";

                    std::cerr << " -- found: " << hn_found << " fn unused: " << fn_unused << " hns: " << het_nodes.size() 
                              << " hom_conflict: " << hom_conflict << "\n";
#endif
                    // covering all hets?
                    if(hn_found + fn_unused == het_nodes.size() && !hom_conflict)
                    {
                        DiploidRef r = {
                            p1 != p2, 
                            false,
                            target[p1].seq(start, end), 
                            target[p2].seq(start, end),
                            refsq
                        };

                        /** homref in the het case means that we have a het+homref call vs. 
                            a het+het call. the two sequences are different by design of the 
                            GraphReference::enumeratePaths function, so if one of them is 
                            equal to the reference, the other one can't be.
                         */
                        if(r.h1 == refsq || r.h2 == refsq)
                        {
                            r.homref = true;
                        }
#ifdef _DEBUG_DIPLOIDREFERENCE
                        std::cerr << "Found pair: " << r.h1 << ":" << r.h2 << " - homref: " << r.homref << " het : " << r.het << "\n";
#endif
                        _impl->di_haps.push_back(r);
                    }
                }
            }
        }
        else
        {
            // no het nodes? all hom -> technically, this should only
            // give us one result, except if filtered variant calls give multiple hom alts
            for(Haplotype & hap : target)
            {
                DiploidRef r = {false, hap.noVar(), hap.seq(start, end), "", refsq};
                if(r.h1 == refsq)
                {
                    r.homref = true;
                }
                _impl->di_haps.push_back(r);
            }
        }
    }
    _impl->pos = _impl->di_haps.begin();
}

/**
 * Get next pair of haplotype sequences
 */
bool DiploidReference::hasNext()
{
    return _impl->pos != _impl->di_haps.end();
}

DiploidRef & DiploidReference::next()
{
    return *_impl->pos;
}

void DiploidReference::advance()
{
    ++_impl->pos;
}

std::list<DiploidRef> const & DiploidReference::result()
{
    return _impl->di_haps;
}

}  // namespace haplotypes
