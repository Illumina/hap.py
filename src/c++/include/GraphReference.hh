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
 * \brief Graph reference class to store diploid haplotype alternatives
 *
 * \file GraphReference.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Haplotype.hh"
#include "Fasta.hh"
#include "Variant.hh"

#include <vector>

namespace haplotypes
{

struct ReferenceNode;
struct ReferenceEdge;

struct GraphReferenceImpl;
class GraphReference
{
public:
    GraphReference(
        const char * vcfname, 
        const char * sample,
        const char * ref_fasta);
    
    GraphReference(GraphReference const &);
    GraphReference const & operator=(GraphReference const &);
    ~GraphReference();

    std::string const & getVCFName() const;
    std::string const & getSample() const;
    std::string const & getFastaName() const;

    /**
     * @brief Apply filters in the input VCF
     * 
     */
    bool getApplyFilters() const;
    void setApplyFilters(bool filters=false);

    /**
     * Get reference fasta file
     */
    FastaFile & getRefFasta();
    FastaFile const & getRefFasta() const;

    /**
     * @brief Enumerate all possible paths, return haplotype candidates
     * 
     */
    void enumeratePaths(
        std::vector<Haplotype> & target,
        const char * chr,
        int64_t start,
        int64_t end,
        int max_n_paths=-1
    );

    /**
     * Create a diploid reference graph for a region of a VCF
     * Optionally, will return the number of unphased het variants
     * (proxy for number of paths)
     */
    void makeGraph(
        const char * chr,
        int64_t start,
        int64_t end,
        std::vector<ReferenceNode> & nodes,
        std::vector<ReferenceEdge> & edges,
        size_t * nhets=NULL
    );

    /**
     * Create a diploid reference graph from a list of Variants records
     * Optionally, will return the number of unphased het variants
     * (proxy for number of paths)
     */
    void makeGraph(
        std::list<variant::Variants> const & input, int ix,
        std::vector<ReferenceNode> & nodes,
        std::vector<ReferenceEdge> & edges,
        size_t * nhets=NULL
    );

    /**
     * Enumerate paths for a graph
     * 
     * Core function for traversing the reference graph
     * 
     * Parameters:
     * 
     * chr, start and end are used to restrict the range of traversal and to make haplotypes unique by their representation
     * in that range.
     * 
     * nodes and edges are the graph as created by makeGraph()
     * 
     * target receives the enumerated haplotype blocks
     * 
     * Graph traversal can be controlled by specifying a source and a sink (this will only enumerate paths starting at source
     * and optionally ending at sink), and by limiting the number of paths that are output.
     * 
     * Finally, each Haplotype block corresponds to traversing a set of nodes, which can be returned for each haplotype 
     * in the nodes_used vector. This is used in DiploidReference to find the pairs of haplotype blocks in a diploid setting
     * (pairs of haplotypes that cover the whole graph).
     * 
     */
    void enumeratePaths(
        const char * chr,
        int64_t start,
        int64_t end,
        std::vector<ReferenceNode> const & nodes,
        std::vector<ReferenceEdge> const & edges,
        std::vector<Haplotype> & target,
        size_t source=0,
        int max_n_paths=-1,
        size_t sink=(size_t)-1,
        std::vector<std::string> * nodes_used = NULL
    );

private:
    GraphReferenceImpl * _impl;
};


struct ReferenceNode : public variant::RefVar
{
    typedef enum _nodetype_t
    {
        source,
        sink,
        homref,
        alternative,
        invalid
    } nodetype_t;

    typedef enum _color_t
    {
        grey = 0,
        black = 1,
        red = 2,
        blue = 3,
    } color_t;

    ReferenceNode(const char * _chr=NULL, int64_t _start = -1, int64_t _end = -1, 
                  ReferenceNode::nodetype_t _type = ReferenceNode::invalid,
                  ReferenceNode::color_t _color = ReferenceNode::grey) :
        chr(_chr == NULL ? "unknown" : _chr), filtered(false), het(false), 
        type(_type), color(_color)
    {
        start = _start;
        end = _end;
    }

    /**
     * @brief Append a ReferenceNode to a haplotype block
     * 
     * @param ht Haplotype block with end() < start
     */
    void appendToHaplotype(Haplotype & ht)  const;

    std::string chr;

    /* we keep track of filtered calls since these may 
     * overlap with other calls */
    bool filtered;

    /* remember if nodes are half of a het or hetalt call  */
    bool het;
    
    nodetype_t type;
    color_t color;
};

struct ReferenceEdge
{   
    ReferenceEdge(size_t _u, size_t _v) :
        u(_u), v(_v) {}
    size_t u;
    size_t v;
};

} // namespace haplotypes
