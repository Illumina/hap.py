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
 * \brief Alignment interface
 *
 * \file Alignment.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

#include "RefVar.hh"

struct AlignmentParameters
{
    AlignmentParameters()
    {
        // default: BWA-MEM scores, see
        // http://bio-bwa.sourceforge.net/bwa.shtml
        //  a   c   g   t   n
        const int8_t default_substitution_scores[] = {
            1, -4, -4, -4,  0,  // a
           -4,  1, -4, -4,  0,  // c
           -4, -4,  1, -4,  0,  // g
           -4, -4, -4,  1,  0,  // t
            0,  0,  0,  0,  0,  // n
        };

        memcpy(subs_mat, default_substitution_scores, 25*sizeof(int8_t));
        gapo = 6;
        gape = 1;
    }

    int8_t subs_mat[25];

    uint8_t gapo;
    uint8_t gape;

    int8_t maxScore() const
    {
        int8_t max_score = std::numeric_limits<int8_t>::min();

        for (int i = 0; i < 25; ++i)
        {
            max_score = std::max(max_score, subs_mat[i]);
        }
        return max_score;
    }
};

/**
 * @brief Sequence alignment class interface
 */
class Alignment
{
public:
    virtual ~Alignment();

    virtual void setParameters(AlignmentParameters const & ap) = 0;
    virtual void getParameters(AlignmentParameters & ap) = 0;

    /* get the best possible score for comparing two sequences of length len */
    virtual int bestScore(int len);

    /**
     * @brief set target sequence
     */
    virtual void setRef(const char * seq) = 0;
    /**
     * @brief set query sequence
     */
    virtual void setQuery(const char * seq) = 0;
    /**
     * @brief Get the alignment score
     */
    virtual int getScore() = 0;
    /**
     * @brief Get a human-readable cigar string + start + end
     */
    virtual void getCigar(int &r0, int & r1, int &a0, int & a1, std::string &cig) = 0;

    /**
     * @brief Get int* cigar string + start + end
     */
    virtual void getCigar(int &r0, int & r1, int &a0, int & a1, int & n_cigar, uint32_t *& cigar) = 0;

    /**
     * @brief Debug dump, optional
     */
    virtual void dump()
    {
        int r0, r1, a0, a1;
        std::string cig;
        getCigar(r0, r1, a0, a1, cig);
        std::cerr
            << "Score: " << getScore() << "\n"
            << "ref: " << r0 << "-" << r1 << "\n"
            << "alt: " << a0 << "-" << a1 << "\n"
            << "Cigar: " << cig << "\n";
    }

};

/** Alignment result storage */
struct AlignmentResult
{
    AlignmentResult(Alignment * aln)
    {
        score = aln->getScore();
        uint32_t * tmp;
        aln->getCigar(s1, e1, s2, e2, n_cigar, tmp);
        cigar = new uint32_t[n_cigar];
        memcpy(cigar, tmp, sizeof(uint32_t)*n_cigar);
    }

    AlignmentResult() : cigar(NULL) {}
    ~AlignmentResult() { if(cigar) { delete [] cigar; } }
    AlignmentResult(AlignmentResult const & rhs) :
        score(rhs.score), hap1(rhs.hap1), hap2(rhs.hap2),
        s1(rhs.s1), e1(rhs.e1), s2(rhs.s2), e2(rhs.e2),
        n_cigar(rhs.n_cigar)
    {
        if(rhs.n_cigar > 0)
        {
            cigar = new uint32_t[rhs.n_cigar];
            memcpy(cigar, rhs.cigar, sizeof(uint32_t)*rhs.n_cigar);
        }
        else
        {
            cigar = NULL;
        }
    }

    AlignmentResult const & operator=(AlignmentResult const & rhs)
    {
        if(&rhs == this)
        {
            return *this;
        }
        if(cigar)
        {
            delete [] cigar;
        }
        score = rhs.score;
        hap1 = rhs.hap1;
        hap2 = rhs.hap2;
        s1 = rhs.s1;
        e1 = rhs.e1;
        s2 = rhs.s2;
        e2 = rhs.e2;
        n_cigar = rhs.n_cigar;
        if(rhs.n_cigar > 0)
        {
            cigar = new uint32_t[rhs.n_cigar];
            memcpy(cigar, rhs.cigar, sizeof(uint32_t)*rhs.n_cigar);
        }
        else
        {
            cigar = NULL;
        }
        return *this;
    }

    int score;
    std::string hap1;
    std::string hap2;
    int s1, e1, s2, e2;
    uint32_t * cigar;
    int n_cigar;
};



/**
 * @brief Factory interface. Caller has to delete the returned
 * pointer
 */
extern Alignment * makeAlignment(const char * type);

/**
 * @brief Format int encoded Cigar string
 *
 * @param tb for padding with "S" : begin
 * @param te for padding with "S" : end
 * @param altlen for padding with "S" : length of alternate sequence
 * @param n_cigar length of cigar
 * @param cigar int* to cigar entries
 *
 */
std::string makeCigar(int tb, int te, int altlen, int n_cigar, uint32_t * cigar);

/** make variants from a cigar string */
extern void getVariantsFromCigar(std::string const & ref, std::string const & alt,
                                 int r0, int a0,
                                 uint32_t * cigar, int n_cigar,
                                 std::list<variant::RefVar> & target);

/** get stats from a cigar string */
extern void getCigarStats(std::string const & ref, std::string const & alt,
                          int r0, int a0,
                          uint32_t * cigar, int n_cigar,
                          int & softclipped, int & matches,
                          int & mismatches, int & ins, int & del);


/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param aln the alignment interface to use
 * @param vars the primitive records
 */
void realignRefVar(FastaFile const & f, const char * chr, variant::RefVar const & rv, Alignment * aln,
                   std::list<variant::RefVar> & vars);

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param aln the alignment interface to use
 * @param snps the number of snps
 * @param ins the number of insertions
 * @param dels the number of deletions
 * @param homref the number of calls with no variation
 */
void realignRefVar(FastaFile const & f, const char * chr, variant::RefVar const & rv, Alignment * aln,
                   size_t & snps, size_t & ins, size_t & dels, size_t & homref,
                   size_t& transitions, size_t& transversions);


