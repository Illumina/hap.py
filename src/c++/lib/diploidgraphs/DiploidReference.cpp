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
#include <bitset>
#include <map>
#include <unordered_map>

#include "Error.hh"

/* #define _DEBUG_DIPLOIDREFERENCE */

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
    DiploidReferenceImpl(GraphReference const & _gr) :
        chr(""), start(-1), end(-1), max_n_paths(-1),
        gr(_gr)
        {}

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
    _impl(new DiploidReferenceImpl(gr))
{
}

DiploidReference::DiploidReference(DiploidReference const & rhs) : _impl(new DiploidReferenceImpl(rhs._impl->gr))
{
}

DiploidReference const & DiploidReference::operator=(DiploidReference const & rhs)
{
    if(&rhs == this)
    {
        return *this;
    }
    delete _impl;

    DiploidReferenceImpl * t_impl = new DiploidReferenceImpl(rhs._impl->gr);
    _impl = t_impl;

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
    std::list<variant::Variants> const & vars,
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
        _impl->gr.makeGraph(vars, sample_ix, nodes, edges, &nhets);

#ifdef _DEBUG_DIPLOIDREFERENCE
        std::cerr << "nhets : " << nhets << " mp: " << _impl->max_n_paths << "\n";
#endif

        if(_impl->max_n_paths > 0 && pow(2.0, (double)nhets) > _impl->max_n_paths)
        {
            error("Too many het nodes (%i) at %s:%i-%i", nhets, chr, start, end);
        }

        std::string refsq = _impl->gr.getRefFasta().query(chr, start, end);

        std::vector<Haplotype> target;
        std::vector<uint64_t> nodes_used;

        _impl->gr.enumeratePaths(chr, start, end,
                          nodes, edges, target,
                          0, (size_t)-1, _impl->max_n_paths,
                          &nodes_used, &nhets);

#ifdef _DEBUG_DIPLOIDREFERENCE
        std::cerr << "Nodes: " << "\n";
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            std::cerr << i << ": " << nodes[i] << " het:" << nodes[i].het << " filtered: " << nodes[i].filtered << " type: " << nodes[i].type << "\n";
        }
        std::cerr << "Edges: " << "\n";
        for (size_t i = 0; i < edges.size(); ++i)
        {
            std::cerr << i << ": " << edges[i].u << " -> " << edges[i].v << "\n";
        }
        std::cerr << "Haplotypes seen: " << "\n";
        for (size_t i = 0; i < target.size(); ++i)
        {
            std::cerr << target[i].repr(start, end) << "\n";
            std::cerr << "NU: " << std::bitset<64>(nodes_used[i]) << "\n";
        }
#endif

        // if we have het nodes, each pair of haplotypes must cover them all
        if(nhets != 0)
        {
            std::unordered_map<uint64_t, size_t> nu_haps;
            uint64_t nh_mask = 1;
            nh_mask = (nh_mask << (int)(nhets)) - 1;

            for (size_t p1 = 0; p1 < nodes_used.size(); ++p1)
            {
                nu_haps[nodes_used[p1]] = p1;
            }

            while(!nu_haps.empty())
            {
                size_t p1 = nu_haps.begin()->second;

                auto opposite_path = nu_haps.find((~nodes_used[p1]) & nh_mask);
                if(opposite_path != nu_haps.end() && opposite_path != nu_haps.begin())
                {
                    size_t p2 = opposite_path->second;
                    // make order reproducible since map is not ordered
                    if(p2 > p1)
                    {
                        std::swap(p1, p2);
                    }

                    nu_haps.erase(nu_haps.begin());
                    nu_haps.erase(opposite_path);

                    DiploidRef r = {
                        true,
                        false,
                        target[p1].seq(start, end),
                        target[p2].seq(start, end),
                        refsq
                    };

                    /** homref in the het case means that we have a het+homref call vs.
                        a het+het call. the two sequences are different by design of the
                        GraphReference::enumeratePaths function, so if one of them is
                        equal to the reference, the other one shouldn't be.
                     */
                    if( (r.h1 == refsq || r.h2 == refsq) && r.h1 != r.h2)
                    {
                        r.homref = true;
                    }

#ifdef _DEBUG_DIPLOIDREFERENCE
                    std::cerr << "Found pair: " << r.h1 << ":" << r.h2 << " - homref: " << r.homref << " het : " << r.het << "\n";
#endif
                    _impl->di_haps.push_back(r);
                }
                else
                {
#ifdef _DEBUG_DIPLOIDREFERENCE
                    std::cerr << "het path " << nodes_used[p1] << " does not have corresponding opposite path at " <<
                                 chr << ":" << start << "-" << end << std::endl;
#endif
                    nu_haps.erase(nu_haps.begin());
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
    if(_impl->di_haps.empty())
    {
        error("Cannot find matching haplotype pairs at %s:%i-%i", chr, start, end);
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
