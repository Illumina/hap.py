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
 * \brief Wrapper class that creates a Variant Processor and various processing steps
 *
 * \file VariantInput.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "VariantInput.hh"

#include "variant/VariantAlleleRemover.hh"
#include "variant/VariantAlleleSplitter.hh"
#include "variant/VariantTee.hh"
#include "variant/VariantAlleleNormalizer.hh"
#include "variant/VariantLocationAggregator.hh"
#include "variant/VariantHomrefSplitter.hh"
#include "variant/VariantAlleleUniq.hh"
#include "variant/VariantCallsOnly.hh"
#include "variant/VariantPrimitiveSplitter.hh"
#include "variant/VariantLeftPadding.hh"

#include <memory>

namespace variant {

struct VariantInput::VariantInputImpl {
    std::string ref_fasta;
    VariantProcessor proc;

    std::unique_ptr<VariantHomrefSplitter> p_homref_splitter;
    std::unique_ptr<VariantTee> p_tee;
    std::unique_ptr<VariantTee> p_tee_raw;
    std::unique_ptr<VariantProcessor> p_proc_homref;
    std::unique_ptr<VariantProcessor> p_proc_raw;
    std::shared_ptr<VariantFilterInterface> p_filter;
    std::unique_ptr<VariantInjectorStep> p_homref_filter;
    std::unique_ptr<VariantLocationAggregator> p_homref_merger;
    std::unique_ptr<VariantCallsOnly> p_calls_only;
    std::unique_ptr<VariantAlleleRemover> p_allele_remover;
    std::unique_ptr<VariantAlleleNormalizer> p_allele_normalizer;
    std::unique_ptr<VariantLocationAggregator> p_merger;
    std::unique_ptr<VariantAlleleUniq> p_allele_uniq;
    std::unique_ptr<VariantAlleleSplitter> p_allele_splitter;
    std::unique_ptr<VariantPrimitiveSplitter> p_primitive_splitter;
    std::unique_ptr<VariantLocationAggregator> p_merger_p;
    std::unique_ptr<VariantLeftPadding> p_padding;
};

VariantInput::VariantInput(
    const char * _ref_fasta,
    bool leftshift,
    bool refpadding,
    bool trimalleles,
    bool splitalleles,
    int mergebylocation,
    bool uniqalleles,
    bool calls_only,
    bool homref_split,
    bool primitive_split,
    bool homref_output,
    int64_t leftshift_limit,
    bool collect_raw
    ) : _impl(new VariantInputImpl())
{
    _impl->ref_fasta = _ref_fasta;

    if(collect_raw)
    {
        _impl->p_tee_raw = std::move(std::unique_ptr<VariantTee>(new VariantTee()));
        _impl->proc.addStep(*_impl->p_tee_raw);

        _impl->p_proc_raw = std::move(std::unique_ptr<VariantProcessor>(new VariantProcessor()));
        _impl->p_proc_raw->addStep(_impl->p_tee_raw->secondOutput());
    }

    if (homref_split)
    {
        _impl->p_homref_splitter = std::move(std::unique_ptr<VariantHomrefSplitter>(new VariantHomrefSplitter()));
        _impl->proc.addStep(*_impl->p_homref_splitter);

        if(homref_output)
        {
            _impl->p_tee = std::move(std::unique_ptr<VariantTee>(new VariantTee()));
            _impl->proc.addStep(*_impl->p_tee);

            _impl->p_proc_homref = std::move(std::unique_ptr<VariantProcessor>(new VariantProcessor()));
            _impl->p_proc_homref->addStep(_impl->p_tee->secondOutput());

            _impl->p_homref_filter = std::move(std::unique_ptr<VariantInjectorStep>(new VariantInjectorStep()));

            struct hr_filter : public VariantFilterInterface
            {
                bool test(Variants const & v) { return v.anyHomref(); }
            };

            _impl->p_filter = std::move(std::shared_ptr<VariantFilterInterface>(new hr_filter()));
            _impl->p_homref_filter->addFilter(_impl->p_filter);
            _impl->p_proc_homref->addStep(*_impl->p_homref_filter);

            _impl->p_homref_merger = std::move(std::unique_ptr<VariantLocationAggregator>(new VariantLocationAggregator()));
            _impl->p_proc_homref->addStep(*_impl->p_homref_merger);
        }
    }

    if (trimalleles)
    {
        _impl->p_allele_remover = std::move(std::unique_ptr<VariantAlleleRemover>(new VariantAlleleRemover()));
        _impl->proc.addStep(*_impl->p_allele_remover);
    }

    if (splitalleles)
    {
        _impl->p_allele_splitter = std::move(std::unique_ptr<VariantAlleleSplitter>(new VariantAlleleSplitter()));
        _impl->proc.addStep(*_impl->p_allele_splitter);
    }

    if (leftshift && !primitive_split)
    {
        _impl->p_allele_normalizer = std::move(std::unique_ptr<VariantAlleleNormalizer>(new VariantAlleleNormalizer()));
        _impl->p_allele_normalizer->setReference(_ref_fasta);
        _impl->p_allele_normalizer->setEnableRefPadding(refpadding);
        _impl->p_allele_normalizer->setLeftshiftLimit(leftshift_limit);
        _impl->p_allele_normalizer->setEnableHomrefVariants(true);
        _impl->proc.addStep(*_impl->p_allele_normalizer);
    }

    if (mergebylocation && !primitive_split)
    {
        _impl->p_merger = std::move(std::unique_ptr<VariantLocationAggregator>(new VariantLocationAggregator()));
        if (mergebylocation == 2)
        {
            _impl->p_merger->setAggregationType(VariantLocationAggregator::aggregate_hetalt);
        }
        else if (mergebylocation == 3)
        {
            _impl->p_merger->setAggregationType(VariantLocationAggregator::aggregate_ambigous);
        }
        _impl->proc.addStep(*_impl->p_merger);
    }

    if(primitive_split)
    {
        _impl->p_primitive_splitter = std::move(std::unique_ptr<VariantPrimitiveSplitter>(new VariantPrimitiveSplitter()));
        _impl->p_primitive_splitter->setReference(_ref_fasta);
        _impl->proc.addStep(*_impl->p_primitive_splitter);

        if (leftshift)
        {
            _impl->p_allele_normalizer = std::move(std::unique_ptr<VariantAlleleNormalizer>(new VariantAlleleNormalizer()));
            _impl->p_allele_normalizer->setReference(_ref_fasta);
            _impl->p_allele_normalizer->setEnableRefPadding(refpadding);
            _impl->p_allele_normalizer->setLeftshiftLimit(leftshift_limit);
            _impl->p_allele_normalizer->setEnableHomrefVariants(true);
            _impl->proc.addStep(*_impl->p_allele_normalizer);
        }

        _impl->p_merger_p = std::move(std::unique_ptr<VariantLocationAggregator>(new VariantLocationAggregator()));
        _impl->p_merger_p->setAggregationType(VariantLocationAggregator::aggregate_hetalt);
        _impl->proc.addStep(*_impl->p_merger_p);
    }

    if (uniqalleles)
    {
        _impl->p_allele_uniq = std::move(std::unique_ptr<VariantAlleleUniq>(new VariantAlleleUniq()));
        _impl->proc.addStep(*_impl->p_allele_uniq);
    }

    if(refpadding)
    {
        _impl->p_padding = std::move(std::unique_ptr<VariantLeftPadding>(new VariantLeftPadding()));
        _impl->p_padding->setReference(_ref_fasta);
        _impl->proc.addStep(*_impl->p_padding);
    }

    if (calls_only)
    {
        _impl->p_calls_only = std::move(std::unique_ptr<VariantCallsOnly>(new VariantCallsOnly()));
        _impl->proc.addStep(*_impl->p_calls_only);
    }
}

VariantInput::~VariantInput()
{
    delete _impl;
}

/** Read variants from a block into a list */
void VariantInput::get(const char * _chr, int64_t start, int64_t end, std::list<Variants> & vlist)
{
    std::string chr;
    if(!_chr)
    {
        chr = "";
    }
    else
    {
        chr = _chr;
    }

    _impl->proc.rewind(chr.c_str(), start);
    if (_impl->p_tee)
    {
        _impl->p_tee->flush();
    }

    while(_impl->proc.advance())
    {
        Variants & vars = _impl->proc.current();
        if(chr == "")
        {
            chr = vars.chr;
        }
        else if(chr != vars.chr)
        {
            // break on change of chr
            break;
        }
        if(end >= 0 && vars.pos > end)
        {
            break;
        }

        vlist.push_back(vars);
    }
}

/** Direct access to variant processor */
VariantProcessor & VariantInput::getProcessor(processor_id id)
{
    switch(id)
    {
        case raw:
            return *(_impl->p_proc_raw);
        case homref:
            return *(_impl->p_proc_homref);
        default:
            return _impl->proc;
    }
}


}
