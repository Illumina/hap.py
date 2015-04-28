// -*- mode: c++; indent-tabs-mode: nil; -*-
//
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
 *  \brief Helpers for graphs
 *
 * \file GraphUtil.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <stdexcept>

namespace graphutil
{

/** default formatter does nothing */
template <typename _t> 
struct dotFormatter
{
    void operator()(std::ostream & , _t const & )
    {
    }
};

template <typename _e> 
struct dotEdgeTyper
{
    std::string operator()(_e const & )
    {
        return " -> ";
    }
};


/**
 * Generic dot writing. 
 * 
 * Edge class _e must have u and v elements of type size_t which 
 * point into nodes vector.
 * 
 */
template<typename _n, typename _e, 
         typename _node_formatter=dotFormatter<_n>, 
         typename _edge_formatter=dotFormatter<_e>,
         typename _edge_typer=dotEdgeTyper<_e>
        >
static inline void writeGraphDot(
    std::vector<_n> const & nodes,
    std::vector<_e> const & edges,
    const char * filename
)
{
    _node_formatter nf;
    _edge_formatter ef;
    _edge_typer et;
    
    std::ofstream output (filename);

    output << "digraph G { \n";

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        output << "node_" << i << " " << " [";
        nf(output, nodes[i]);
        output << "] ;\n";
    }

    output << "\n";

    for (size_t i = 0; i < edges.size(); ++i)
    {
        output << "node_" << edges[i].u << et(edges[i])
               << "node_" << edges[i].v << " [";
        ef(output, edges[i]);
        output << "] ;\n";
    }
    output << "}\n";
}

/**
 * Make an adjacency list from node/edge list
 */
template<typename _e>
void adjList(
    size_t nnodes,
    std::vector<_e> const & edges,
    std::vector< std::list< size_t > > & adj
    )
{
    adj.resize(nnodes);
    for (size_t i = 0; i < adj.size(); ++i)
    {
        adj[i].clear();
    }

    for (_e const & edge : edges)
    {
        adj[edge.u].push_back(edge.v);
    }
}

/**
 * Make an inverse adjacency list from node/edge list (reverse edges)
 */
template<typename _e>
void invAdjList(
    size_t nnodes,
    std::vector<_e> const & edges,
    std::vector< std::list< size_t > > & adj
    )
{
    adj.resize(nnodes);
    for (size_t i = 0; i < adj.size(); ++i)
    {
        adj[i].clear();
    }

    for (_e const & edge : edges)
    {
        adj[edge.v].push_back(edge.u);
    }
}

/**
 * Compute in-degrees for all nodes
 */
template<typename _e>
void inDegrees(
    size_t nnodes,
    std::vector<_e> const & edges,
    std::vector< size_t > & indegrees
    )
{
    size_t zero = 0;
    indegrees.resize(nnodes, zero);

    for (_e const & edge : edges)
    {
        indegrees[edge.v]++;
    }
}

/**
 * Compute in-degrees for all nodes
 */
template<typename _e>
void outDegrees(
    size_t nnodes,
    std::vector<_e> const & edges,
    std::vector< size_t > & outdegrees
    )
{
    size_t zero = 0;
    outdegrees.resize(nnodes, zero);

    for (_e const & edge : edges)
    {
        outdegrees[edge.u]++;
    }
}

}   // namespace graphutil

