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
 * \brief Graph reference implementation
 *
 * \file GraphReference.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "GraphReference.hh"
#include "Variant.hh"
#include "Haplotype.hh"
#include "Fasta.hh"
#include "helpers/StringUtil.hh"
#include "helpers/GraphUtil.hh"

#include "variant/VariantAlleleRemover.hh"
#include "variant/VariantAlleleSplitter.hh"
#include "variant/VariantAlleleNormalizer.hh"
#include "variant/VariantLocationAggregator.hh"
#include "variant/VariantAlleleUniq.hh"

#include <cassert>
#include <memory>
#include <list>
#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <bitset>

// #define _DEBUG_GRAPHREFERENCE

#ifndef MAX_GRAPHREFERENCE_NODES
#define MAX_GRAPHREFERENCE_NODES 4096
#endif

#include "Error.hh"


using namespace variant;

namespace haplotypes
{

struct GraphReferenceImpl
{
    GraphReferenceImpl( const char * _ref_fasta) : refsq(_ref_fasta)
    {
    }

    GraphReferenceImpl const & operator=(GraphReferenceImpl const & rhs)
    {
        if(this == &rhs)
        {
            return *this;
        }
        refsq = FastaFile(rhs.refsq.getFilename().c_str());
        return *this;
    }

    FastaFile refsq;
};

GraphReference::GraphReference(const char * ref_fasta) :
        _impl(new GraphReferenceImpl(ref_fasta)) {}

GraphReference::GraphReference(GraphReference const & rhs) :
    _impl(new GraphReferenceImpl(rhs._impl->refsq.getFilename().c_str()))
{
}

GraphReference & GraphReference::operator=(GraphReference const & rhs)
{
    if(this == &rhs)
    {
        return *this;
    }
    *_impl = *rhs._impl;
    return *this;
}

GraphReference::~GraphReference()
{
    delete _impl;
}

FastaFile & GraphReference::getRefFasta()
{
    return _impl->refsq;
}

FastaFile const & GraphReference::getRefFasta() const
{
    return _impl->refsq;
}

/**
 * Create a diploid reference graph from a list of Variants records
 * Optionally, will return the number of unphased het variants
 * (proxy for number of paths)
 */
void GraphReference::makeGraph(
    std::list<variant::Variants> const & input, int ix,
    std::vector<ReferenceNode> & nodes,
    std::vector<ReferenceEdge> & edges,
    size_t * nhets
)
{
    if(nhets != NULL)
    {
        *nhets = 0;
    }
    nodes.clear();
    edges.clear();

#ifdef _DEBUG_GRAPHREFERENCE
    std::cerr << "Building graph from variant set, " << input.size() << " variants / " << " sample-index " << ix << "\n";
#endif

    std::string chr;
    int64_t start = 0;

    if (!input.empty())
    {
        chr = input.front().chr;
        start = input.front().pos;
    }

    std::list<size_t> previous;

    // insert source
    ReferenceNode rn(
        chr.c_str(), start-1, start-1,
        ReferenceNode::source
    );
    nodes.push_back(rn);
    previous.push_back(0);

    ReferenceNode current[2];
    for(Variants const & vars : input)
    {
#ifdef _DEBUG_GRAPHREFERENCE
        std::cerr << stringutil::formatPos(chr, vars.pos) << ": " << vars << "\n";
#endif

        Call const & call = vars.calls[ix];

        for(int j = 0; j < 2; ++j)
        {
            current[j].type = ReferenceNode::invalid;
            current[j].color = ReferenceNode::grey;
            current[j].chr = chr;
            current[j].filtered = false;
            current[j].het = false;
        }

        std::string refal = "N";
        if(vars.variation.size() > 0)
        {
            refal = _impl->refsq.query(chr.c_str(), vars.pos, vars.pos + vars.len - 1);
        }

        switch(getGTType(call))
        {
            case gt_het:
            case gt_hetalt:
                current[0].het = true;
                current[1].het = true;
                if(call.nfilter > 0 && !(call.nfilter == 1 && call.filter[0] == "PASS"))
                {
                    current[0].filtered = true;
                    current[1].filtered = true;
                }
                if(call.phased)
                {
                    current[0].color = ReferenceNode::red;
                    current[1].color = ReferenceNode::blue;
                }
                else if(nhets != NULL)
                {
                    ++*nhets;
                }
                for(int j = 0; j < 2; ++j)
                {
                    if(call.gt[j] == 0)
                    {
                        current[j].type = ReferenceNode::homref;
                        current[j].start = vars.pos;
                        current[j].end = vars.pos + vars.len - 1;
                        current[j].alt = refal;
                    }
                    else
                    {
                        current[j].type = ReferenceNode::alternative;
                        current[j].start = vars.variation[call.gt[j]-1].start;
                        current[j].end = vars.variation[call.gt[j]-1].end;
                        current[j].alt = vars.variation[call.gt[j]-1].alt;
                    }
                }
                break;
            case gt_homref:
                // do not generate homref nodes -- when they overlap, things go bad downstream
                // also, we don't really use homref information in the graph.
                // current[0].color = ReferenceNode::black;
                // current[0].type = ReferenceNode::homref;
                // current[0].start = vars.pos;
                // current[0].end = vars.pos + vars.len - 1;
                // current[0].alt = refal;
                break;
            case gt_haploid:
                // TODO: we might want to pass on some information into the reference graph here
            case gt_homalt:
                if(call.nfilter > 0 && !(call.nfilter == 1 && call.filter[0] == "PASS"))
                {
                    current[0].filtered = true;
                }
                current[0].color = ReferenceNode::black;
                current[0].type = ReferenceNode::alternative;
                current[0].start = vars.variation[call.gt[0]-1].start;
                current[0].end = vars.variation[call.gt[0]-1].end;
                current[0].alt = vars.variation[call.gt[0]-1].alt;
                current[0].het = false;
                break;
            default:
                // all other records are ignored
                break;
        }

        // TODO: We can probably do better than reassigning each time.
        std::list<size_t> next_previous;
        size_t current_pos[2];
        for (int j = 0; j < 2; ++j)
        {
            if(current[j].type != ReferenceNode::invalid)
            {
                if(nodes.size() >= MAX_GRAPHREFERENCE_NODES)
                {
                    error("Cannot create reference graph with more than %i nodes (this is determined at compile time).", MAX_GRAPHREFERENCE_NODES);
                }
                current_pos[j] = nodes.size();
                nodes.push_back(current[j]);
                next_previous.push_back(current_pos[j]);
            }
            else
            {
                current_pos[j] = (size_t)-1;
            }
        }

        bool has_predecessor[2] = {false,false};

        // connect the nodes
        for(size_t x : previous)
        {
            ReferenceNode & node(nodes[x]);

            // will be set to true if we found at least one new node to follow this one
            bool superseeded = false;

            for(int j = 0; j < 2; ++j)
            {
                if(current_pos[j] == (size_t) -1)
                {
                    continue;
                }

                // to connect, the previous node must end before the current one starts, and
                // the colors must be "compatible":
                // Black and grey (both haplotypes, known/unknown phasing) connect to everything
                // Red and blue do not connect to each other
                if( node.type == ReferenceNode::source ||
                   (node.end < current[j].start &&
                   (   node.color <= ReferenceNode::black
                    || current[j].color <= ReferenceNode::black
                    || node.color == current[j].color) ) )
                {
                    has_predecessor[j] = true;
                    edges.push_back(ReferenceEdge(x, current_pos[j]));
                    superseeded = true;
                }
            }

            if(!superseeded)
            {
                next_previous.push_back(x);
            }
        }
        previous = next_previous;

        // if variants are out of order, we might need to
        // find a predecessor
        for(int j = 0; j < 2; ++j)
        {
            if(current_pos[j] != (size_t) -1 && !has_predecessor[j])
            {
                ReferenceNode & node = nodes[current_pos[j]];
                // search backwards
                size_t k = current_pos[j]-1;

                while(k != (size_t) -1)
                {
                    if(nodes[k].end < node.start &&
                       (   node.color <= ReferenceNode::black
                        || nodes[k].color <= ReferenceNode::black
                        || node.color == nodes[k].color)
                    )
                    {
                        edges.push_back(ReferenceEdge(k, current_pos[j]));
                        has_predecessor[j] = true;
                        break;
                    }
                    --k;
                }
                // technically, we should have gotten here in the loop above already.
                if(!has_predecessor[j])
                {
                    edges.push_back(ReferenceEdge(0, current_pos[j]));
                }
            }
        }
    }

    // insert sink
    int64_t max_end = 0;
    for(size_t x : previous)
    {
        ReferenceNode & node(nodes[x]);
        max_end = std::max(node.end, max_end);
        edges.push_back(ReferenceEdge(x, nodes.size()));
    }
    ReferenceNode rs(
        chr.c_str(), max_end, max_end,
        ReferenceNode::sink
    );
    nodes.push_back(rs);
}


/**
 * Enumerate paths for a reference graph
 */
void GraphReference::enumeratePaths(
    const char * chr,
    int64_t start,
    int64_t end,
    std::vector<ReferenceNode> const & nodes,
    std::vector<ReferenceEdge> const & edges,
    std::vector<Haplotype> & target,
    size_t source,
    int max_n_paths,
    size_t sink,
    std::vector<std::string> * nodes_used_vec
)
{
    // make adjacency list from edge list
    std::vector< std::list< size_t > > adj;
    graphutil::adjList(nodes.size(), edges, adj);

#ifdef _DEBUG_GRAPHREFERENCE
    static int bp_id_ctr = 0;
#endif

    typedef std::map<std::string, std::unique_ptr<Haplotype> > unique_hap_map_t;
    typedef std::map<std::string, std::string > hapinfo_map_t;

    // we assume the reference graph is loop free
    // otherwise, this doesn't really work
    typedef struct _branchpoint
    {
        _branchpoint(size_t _node,
                     std::list<size_t>::iterator _next_choice,
                     unique_hap_map_t::iterator _up_to_here,
                     std::bitset<MAX_GRAPHREFERENCE_NODES> const & _nodes_used) :
            node(_node), color(ReferenceNode::grey), next_choice(_next_choice),
            up_to_here(_up_to_here), nodes_used(_nodes_used)
#ifdef _DEBUG_GRAPHREFERENCE
            , bp_id(bp_id_ctr++)
#endif
        {}

        size_t node;
        ReferenceNode::color_t color;
        std::list<size_t>::iterator next_choice;
        unique_hap_map_t::iterator up_to_here; // pointer into unique haplotype block

        // track the number of nodes we used
        std::bitset<MAX_GRAPHREFERENCE_NODES> nodes_used;

#ifdef _DEBUG_GRAPHREFERENCE
        int bp_id;
#endif
    } branchpoint;

    std::list< branchpoint > hlist;

    // we contract the variation graph into a "unique haplotype block graph"
    // this map gives the nodes (indexed by canonical haplotype representation).
    unique_hap_map_t unique_haps;
    hapinfo_map_t hapinfos;

    // list of hap ids we get after finishing traversal
    std::set<std::string> final_haps;

    // generate starting hap block
    std::unique_ptr<Haplotype> ph0 (new Haplotype(chr, _impl->refsq.getFilename().c_str()));
    Haplotype & h0(*ph0);
    nodes[source].appendToHaplotype(h0);
    std::string current_rp(h0.repr(start, end));
    std::string first_rp(current_rp);

    std::bitset<MAX_GRAPHREFERENCE_NODES> initial_bitset;
    initial_bitset[0] = 1;

    unique_haps.emplace(current_rp, std::move(ph0));
    hapinfos.emplace(current_rp, initial_bitset.to_string().substr(MAX_GRAPHREFERENCE_NODES-nodes.size(), nodes.size()));

    // first branch point
    branchpoint bp(source, adj[source].begin(), unique_haps.find(current_rp), initial_bitset);
    hlist.push_back(bp);

    while(!hlist.empty() && final_haps.size() < ((size_t)max_n_paths))
    {
        branchpoint & current(hlist.front());

        // exhaustively go through all the choices from this point on
        while(current.next_choice != adj[current.node].end())
        {
            size_t nextone = *current.next_choice;
            // remember for next time we come here that we've gone
            // through the first branch
            current.next_choice++;

            // copy unique haplotype to start from
            std::unique_ptr<Haplotype> pht(new Haplotype(*(current.up_to_here->second)));
            Haplotype & ht(*pht);
            std::bitset<MAX_GRAPHREFERENCE_NODES> nodes_used(current.nodes_used);

#ifdef _DEBUG_GRAPHREFERENCE
            std::cerr << "(Re)starting at bp " << current.bp_id << ": "
                                               << current.up_to_here->first
                                               << " / " << current.up_to_here->second->repr() << "\n";
#endif
            bool cont = true;
            ReferenceNode::color_t current_path_color = current.color;
            while(cont)
            {
                if(sink != (size_t)-1 && nodes[nextone].start > nodes[sink].start)
                {
#ifdef _DEBUG_GRAPHREFERENCE
                    std::cerr << "Stopping enumeration because we've passed the sink position at " <<
                        nodes[nextone].start << " / " << nodes[sink].start << "\n";
#endif
                    cont = false;
                    break;
                }

                //    make sure when traversing phased paths we don't switch haplotypes!
                // == make sure paths which contain a red node will never contain a blue one
                //    and vice versa
                //    implemented by checking color compatibility of current path and next node
                if(   nodes[nextone].color > ReferenceNode::black
                   && current_path_color > ReferenceNode::black
                   && nodes[nextone].color != current_path_color)
                {
                    // not compatible => terminate path
#ifdef _DEBUG_GRAPHREFERENCE
                    std::cerr << "Ignoring incompatible branch -- path color is " << current_path_color << " next color is " << nodes[nextone].color << "\n";
#endif
                    cont = false;
                    break;
                }

                // mark that we used this node on this path
                nodes_used[nextone] = 1;

                nodes[nextone].appendToHaplotype(ht);
                // update path color
                current_path_color = std::max(nodes[nextone].color, current_path_color);

#ifdef _DEBUG_GRAPHREFERENCE
                std::cerr << "Appending to HT, now: " << ht.repr() << "\n";
#endif
                current_rp = h0.repr(start, end);

                if(nextone == sink || adj[nextone].size() != 1)
                                         // end of path, Haplotype block is complete
                                         // or more than one choice: insert a haplotype block node if necessary
                                         // and link up to previous
                {
                    std::string ht_rp(ht.repr(start, end));
                    unique_hap_map_t::iterator pos = unique_haps.find(ht_rp);
                    if(pos == unique_haps.end())
                    {
                        // insert a copy
                        std::unique_ptr<Haplotype> pht_copy(new Haplotype(ht));
                        pos = unique_haps.emplace(ht_rp, std::move(pht_copy)).first;
                        hapinfos.emplace(ht_rp, nodes_used.to_string().substr(MAX_GRAPHREFERENCE_NODES-nodes.size(), nodes.size()));
#ifdef _DEBUG_GRAPHREFERENCE
                        std::cerr << "Adding UHB " << ht_rp << "\n";
                    }
                    else
                    {
                        std::cerr << "UHB " << ht_rp << " already exists.\n";
#endif
                    }
                    if(nextone != sink && adj[nextone].size() > 1) // more than 1 choice => create branch point
                    {   // more than one choice: create branchpoint
                        branchpoint bp(nextone, std::next(adj[nextone].begin()), pos, nodes_used);
                        bp.color = current_path_color;
                        hlist.push_back(bp);
#ifdef _DEBUG_GRAPHREFERENCE
                        std::cerr << "Creating BP " << bp.bp_id << "\n";
#endif
                    }
                    else
                    {   // no further options
                        cont = false;
                        if(sink == (size_t)-1 || nextone == sink)
                        {
                            final_haps.insert(ht_rp);
#ifdef _DEBUG_GRAPHREFERENCE
                            std::cerr << "Finished path from BP " << current.bp_id << " at " << ht_rp << "\n";
#endif
                        }
                        else
                        {
#ifdef _DEBUG_GRAPHREFERENCE
                            std::cerr << "Finished path from BP" << current.bp_id << " at " << ht_rp << " (no sink)\n";
#endif
                        }
                        break;
                    }
                }
                nextone = adj[nextone].front();
            }
        }
        hlist.pop_front();
    }

    for(std::string const & p : final_haps)
    {
        // save all blocks with no out edges
        target.push_back(*(unique_haps[p]));
        if(nodes_used_vec != NULL)
        {
            nodes_used_vec->push_back(hapinfos[p]);
        }
    }

    // no final haps because all is homref?
    if(target.empty())
    {
        target.push_back(Haplotype(chr, _impl->refsq.getFilename().c_str()));

        if(nodes_used_vec != NULL)
        {
            nodes_used_vec->push_back("1");
        }
    }
}

/**
 * @brief Append a ReferenceNode to a haplotype block
 *
 * @param ht Haplotype block with end() < start
 */
void ReferenceNode::appendToHaplotype(Haplotype & ht) const
{
    if(type == alternative)
    {
        ht.addVar(start, end, alt);
    }
}

} // namespace haplotypes
