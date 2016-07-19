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
#include <queue>

/* #define _DEBUG_GRAPHREFERENCE */

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
    std::list<variant::Variants> const & _input, int ix,
    std::vector<ReferenceNode> & nodes,
    std::vector<ReferenceEdge> & edges,
    size_t * nhets
)
{
    std::queue<Variants> input;

    // trim and split alleles
    for(auto const & v : _input)
    {
        Call const & c = v.calls[ix];
        switch(getGTType(c))
        {
            case gt_homalt:
            case gt_het:
            {
                Variants v_copy;
                v_copy.id = v.id;
                v_copy.chr = v.chr;
                v_copy.calls.push_back(c);
                int which_al = -1;
                if(v_copy.calls[0].gt[0] != 0)
                {
                    which_al = v_copy.calls[0].gt[0] - 1;
                    v_copy.calls[0].gt[0] = 1;
                }
                if(v_copy.calls[0].gt[1] != 0)
                {
                    which_al = v_copy.calls[0].gt[1] - 1;
                    v_copy.calls[0].gt[1] = 1;
                }
                if(which_al < 0 || which_al >= (int)v.variation.size())
                {
                    error("invalid GT at %s:%i", v.chr.c_str(), v.pos);
                }
                v_copy.variation.push_back(v.variation[which_al]);
                trimLeft(_impl->refsq, v.chr.c_str(), v_copy.variation[0], false);
                trimRight(_impl->refsq, v.chr.c_str(), v_copy.variation[0], false);
                if(v_copy.variation[0].start <= v_copy.variation[0].end)
                {
                    v_copy.pos = v_copy.variation[0].start;
                    v_copy.len = v_copy.variation[0].end - v_copy.variation[0].start + 1;
                }
                else
                {
                    v_copy.pos = v_copy.variation[0].start - 1;
                    v_copy.len = v_copy.variation[0].end - v_copy.variation[0].start + 1;
                }
                input.push(v_copy);
            }
            break;
            case gt_hetalt:
            {
                Variants v_copy;
                v_copy.id = v.id;
                v_copy.chr = v.chr;
                v_copy.calls.resize(1);
                v_copy.calls[0] = c;
                int64_t min_pos = std::numeric_limits<int64_t>::max();
                int64_t max_pos = -1;

                for(size_t g = 0; g < c.ngt; ++g)
                {
                    if(c.gt[g] <= 0 || c.gt[g] > (int)v.variation.size())
                    {
                        error("invalid GT at %s:%i", v.chr.c_str(), v.pos);
                    }
                    v_copy.variation.push_back(v.variation[c.gt[g]-1]);
                    v_copy.calls[0].gt[g] = v_copy.variation.size();
                    trimLeft(_impl->refsq, v.chr.c_str(), v_copy.variation.back(), false);
                    trimRight(_impl->refsq, v.chr.c_str(), v_copy.variation.back(), false);
                    if(min_pos > v_copy.variation.back().start || min_pos > v_copy.variation.back().end)
                    {
                        min_pos = std::min(v_copy.variation.back().start, v_copy.variation.back().end);
                    }
                    if(max_pos < v_copy.variation.back().start || max_pos < v_copy.variation.back().end)
                    {
                        max_pos = std::max(v_copy.variation.back().start, v_copy.variation.back().end);
                    }
                }
                if(min_pos <= max_pos)
                {
                    v_copy.pos = min_pos;
                    v_copy.len = max_pos - min_pos + 1;
                }
                else
                {
                    v_copy.pos = min_pos - 1;
                    v_copy.len = max_pos - min_pos + 1;
                }
                input.push(v_copy);
            }
            break;
            default: break;
        }
    }
    // we have retrieved calls only for our sample
    ix = 0;
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
    // dynamically build adjacency list
    std::vector< std::set <size_t> > adj;
    for(; !input.empty(); input.pop())
    {
        Variants const & vars = input.front();
#ifdef _DEBUG_GRAPHREFERENCE
        std::cerr << "VAR: " << stringutil::formatPos(chr, vars.pos) << ": " << vars << "\n";
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
            case gt_homalt:
                if(call.nfilter > 0 && !(call.nfilter == 1 && call.filter[0] == "PASS"))
                {
                    current[0].filtered = true;
                }
                current[0].color = ReferenceNode::black;
                current[0].type = ReferenceNode::alternative;
                {
                    RefVar rv = vars.variation[call.gt[0]-1];
                    current[0].start = rv.start;
                    current[0].end = rv.end;
                    current[0].alt = rv.alt;
                }
                current[0].het = false;
                break;
            default:
                break;
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
                if(call.gt[0] == 0 || call.gt[1] == 0)
                {
                    int j = 0;
                    if(call.gt[j] == 0)
                    {
                        ++j;
                    }

                    RefVar rv = vars.variation[call.gt[j]-1];

                    current[j].type = ReferenceNode::alternative;
                    current[j].start = rv.start;
                    current[j].end = rv.end;
                    current[j].alt = rv.alt;

                    j = j ^ 1;

                    current[j].type = ReferenceNode::homref;
                    if(rv.start > rv.end)
                    {
                        // insertion => hom-ref block of length 0
                        current[j].start = std::max(rv.start, rv.end);
                    }
                    else
                    {
                        current[j].start = rv.start;
                    }
                    current[j].end = current[j].start - 1;
                    current[j].alt = "";
                }
                else
                {
                    RefVar const & rv1 = vars.variation[call.gt[0]-1];
                    current[0].type = ReferenceNode::alternative;
                    current[0].start = rv1.start;
                    current[0].end = rv1.end;
                    current[0].alt = rv1.alt;

                    RefVar const & rv2 = vars.variation[call.gt[1]-1];
                    current[1].type = ReferenceNode::alternative;
                    current[1].start = rv2.start;
                    current[1].end = rv2.end;
                    current[1].alt = rv2.alt;
                }
                break;
        }

        // List of "previous nodes" to connect
        // TODO: We can probably do better than reassigning each time.
        std::list<size_t> next_previous;

        // determine list ids for the generated nodes
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
        // make sure adjacency list is large enough
        adj.resize(nodes.size());

        bool has_predecessor[2] = {false,false};

        // connect the nodes
        for(size_t x : previous)
        {
            ReferenceNode & node(nodes[x]);

            bool superseeded = false;

            // will be set to true if we found at least one new node to follow this one
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
                    adj[x].insert(current_pos[j]);
                    superseeded = true;
                }
            }

            if(!superseeded)
            {
                next_previous.insert(next_previous.end(), x);
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
                        adj[k].insert(current_pos[j]);
                        has_predecessor[j] = true;
                        break;
                    }
                    --k;
                }
                // technically, we should have gotten here in the loop above already.
                if(!has_predecessor[j])
                {
                    adj[0].insert(current_pos[j]);
                }
            }
        }
    }

    // insert sink
    int64_t max_end = 0;
    adj.resize(nodes.size());
    for(size_t x : previous)
    {
        ReferenceNode & node(nodes[x]);
        max_end = std::max(node.end, max_end);
        adj[x].insert(nodes.size());
    }
    ReferenceNode rs(
        chr.c_str(), max_end, max_end,
        ReferenceNode::sink
    );
    nodes.push_back(rs);

    for(size_t from = 0; from < adj.size(); ++from)
    {
        for(size_t succ : adj[from])
        {
            edges.push_back(ReferenceEdge(from, succ));
        }
    }
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
    size_t sink,
    int max_n_paths,
    std::vector<uint64_t> * nodes_used_vec,
    size_t * n_hets
)
{
    // make adjacency list from edge list
    std::vector< std::list< size_t > > adj;
    graphutil::adjList(nodes.size(), edges, adj);

    std::vector<uint64_t> node_masks(nodes.size());
    size_t hets = 0;
    uint64_t het_mask = 1;
    size_t homs = 0;
    for(size_t ni = 0; ni < nodes.size(); ++ni)
    {
        auto const & n  = nodes[ni];
        node_masks[ni] = 0;

        bool outside_source_sink = ni < source || (sink != (size_t)-1 && n.start > nodes[sink].start);
        if(n.color == ReferenceNode::black && n.type == ReferenceNode::alternative && !outside_source_sink)
        {
            // this is the number of homs we expect between source and sink
            ++homs;
        }
        else if(n.color != ReferenceNode::black && n.type == ReferenceNode::alternative && !outside_source_sink)
        {
            ++hets;
            if(hets > 63)
            {
                error("Maximum number of hets in block exceeded for %s:%i-%i", chr, start, end);
            }
            node_masks[ni] = het_mask;
#ifdef _DEBUG_GRAPHREFERENCE
            std::cerr << "Node " << n << " has mask " << std::bitset<64>(het_mask).to_string() << "\n";
#endif
            het_mask <<= 1;
        }
    }

    if(n_hets)
    {
        *n_hets = hets;
    }

#ifdef _DEBUG_GRAPHREFERENCE
    static int bp_id_ctr = 0;
#endif

    // we assume the reference graph is loop free
    // otherwise, this doesn't really work
    typedef struct _branchpoint
    {
        _branchpoint(size_t _node,
                     std::list<size_t>::iterator _next_choice,
                     ReferenceNode::color_t _color,
                     Haplotype const & _up_to_here,
                     uint64_t _nodes_used,
                     std::set < std::string > const & _sequences_seen,
                     size_t _homs_used
        ) :
            node(_node),  next_choice(_next_choice), color(_color),
            up_to_here(_up_to_here), nodes_used(_nodes_used),
            sequences_seen(_sequences_seen), homs_used(_homs_used)
#ifdef _DEBUG_GRAPHREFERENCE
            , bp_id(bp_id_ctr++)
#endif
        {}

        size_t node;  // node to start with

        std::list<size_t>::iterator next_choice; // next edge to go through

        ReferenceNode::color_t color; // path color for respecting phasing

        Haplotype up_to_here; // observed haplotype up to here

        uint64_t nodes_used; // track which nodes were used

        std::set<std::string> sequences_seen; // HAP-147 track sequences we have seen already

        size_t homs_used;  // count the hom variants we have used already

#ifdef _DEBUG_GRAPHREFERENCE
        int bp_id;
#endif
    } branchpoint;

    std::list< branchpoint > hlist;

    // generate starting hap block
    ReferenceNode::color_t current_path_color = nodes[source].color;
    Haplotype ht(chr, _impl->refsq.getFilename().c_str());
    nodes[source].appendToHaplotype(ht);
    std::set<std::string> sequences_seen;
    sequences_seen.insert(ht.seq(start, end));
    uint64_t nodes_used = node_masks[source];
    size_t homs_used = 0;
    if(nodes[source].color == ReferenceNode::black && nodes[source].type == ReferenceNode::alternative)
    {
        ++homs_used;
    }

    // first branch point
    hlist.push_back(branchpoint(source, adj[source].begin(), current_path_color, ht,
                                nodes_used, sequences_seen, homs_used));

    while(!hlist.empty() && target.size() < ((size_t)max_n_paths))
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
            current_path_color = current.color;
            ht = current.up_to_here;
            nodes_used = current.nodes_used;
            sequences_seen = current.sequences_seen;
            homs_used = current.homs_used;

#ifdef _DEBUG_GRAPHREFERENCE
            std::cerr << "(Re)starting at bp " << current.bp_id << ": "
                                               << current.up_to_here.repr()
                                               << " nu: "
                                               << std::bitset<64>(current.nodes_used).to_string();

            std::cerr << " sequences_seen: ";
            for (auto x : sequences_seen) {
                std::cerr << " " << x;
            }

            std::cerr << "\n";
#endif
            bool cont = true;
            while(cont)
            {

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

                // update path color
                current_path_color = std::max(nodes[nextone].color, current_path_color);

                nodes[nextone].appendToHaplotype(ht);

                // HAP-147: check that we haven't appended a variant that brought us back
                //          to a sequence we already observed (e.g. insert an A and then
                //          delete it again). We assume that VCFs don't contain cases where
                //          this is valid choice of paths.
                // Note we only need to do this if we already used het variants. Otherwise,
                // we don't really have a choice and need to produce the same sequence twice.
                std::string modified_rp = ht.seq(start, end);
                // mark that we used this node on this path
                nodes_used |= node_masks[nextone];
                if(nodes[nextone].type == ReferenceNode::alternative && nodes_used != 0)
                {
                    if(sequences_seen.count(modified_rp))
                    {
#ifdef _DEBUG_GRAPHREFERENCE
                        std::cerr << "Ignoring branch where we see the same sequence twice " << start << "-" << end << ": " << modified_rp << " seen: ";
                        for (auto i : sequences_seen) {
                            std::cerr << i << ", ";
                        }
                        std::cerr << "\n";
#endif
                        cont = false;
                        break;
                    }
                    sequences_seen.insert(modified_rp);
                }

                if(nodes[nextone].color == ReferenceNode::black && nodes[nextone].type == ReferenceNode::alternative)
                {
                    ++homs_used;
                }

#ifdef _DEBUG_GRAPHREFERENCE
                std::cerr << "Appending " << nodes[nextone] << " to HT, now: " << modified_rp << "\n";
#endif

                // end of path?
                if(   nextone == sink || adj[nextone].size() == 0
                   || (sink != (size_t)-1 && nodes[nextone].start > nodes[sink].start))
                {
                    cont = false;
                    if(homs_used == homs && !(sink != (size_t)-1 && nodes[nextone].start > nodes[sink].start))
                    {
                        // save all blocks with no out edges
                        target.push_back(ht);
                        if(nodes_used_vec != NULL)
                        {
                            nodes_used_vec->push_back(nodes_used);
                        }
#ifdef _DEBUG_GRAPHREFERENCE
                        std::cerr << "Finished path from BP " << current.bp_id << " at " << target[target.size()-1].seq(start, end);
                        if(nodes_used_vec)
                        {
                            std::cerr << " u: " << std::bitset<64>((*nodes_used_vec)[nodes_used_vec->size()-1]).to_string();
                        }
                        std::cerr << "\n";
#endif
                    }
                    else
                    {
#ifdef _DEBUG_GRAPHREFERENCE
                        std::cerr << "Finished path from BP" <<
                            current.bp_id << " at " << ht.seq(start, end) <<
                            " ... this path does not cover all expected hom variants / doesn't end at sink"
                            " and will be ignored\n";
#endif
                    }
#ifdef _DEBUG_GRAPHREFERENCE
                    std::cerr << "------------------\n\n";
#endif
                    break;
                }

                // we know that adj[nextone].size() > 0
                if(adj[nextone].size() > 1)
                {   // more than one choice: create branchpoint
                    hlist.push_back(
                        branchpoint(nextone,
                                    std::next(adj[nextone].begin()),
                                    current_path_color,
                                    ht,
                                    nodes_used,
                                    sequences_seen,
                                    homs_used));
#ifdef _DEBUG_GRAPHREFERENCE
                    std::cerr << "Creating BP " << hlist.back().bp_id << "\n";
#endif
                }
                nextone = adj[nextone].front();
            }
        }
        hlist.pop_front();
    }

    // no final haps because all is homref?
    if(target.empty() && hets == 0 && homs == 0)
    {
        target.push_back(Haplotype(chr, _impl->refsq.getFilename().c_str()));

        if(nodes_used_vec != NULL)
        {
            nodes_used_vec->push_back(0);
        }
    }

    if(target.empty())
    {
        error("Failed to create any haplotype sequences from variants at %s:%i-%i", chr, start, end);
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
