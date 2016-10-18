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
 * Variant Processor
 *
 * \file VariantProcessor.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "Variant.hh"

#include <boost/range/adaptor/reversed.hpp>

//#define DEBUG_VARIANTPROCESSOR
//#define DEBUG_VARIANTPROCESSOR_STEPS

namespace variant {

/**
 * @brief Add variants from a VariantReader (see above for VariantBufferMode and param)
 */
bool AbstractVariantProcessingStep::add(VariantReader & source,
        VariantBufferMode bt,
        int64_t param
       )
{
    int64_t end = -1;
    switch(bt)
    {
        case VariantBufferMode::buffer_block:
            end = -2;
            break;
        case VariantBufferMode::buffer_endpos:
            end = param;
            break;
        case VariantBufferMode::buffer_all:
        default:
            break;
    }

    std::string chr = "";
    int64_t last_end = -1;
    int64_t count = 0;
    bool advance_success = false;
    bool any_calls = false;
    while(source.advance() == true)
    {
        Variants & v = source.current();
        // return true if any variants still in source
        advance_success = true;
        if (bt == VariantBufferMode::buffer_count && count+1 > param)
        {
            // put last variant back
            source.enqueue(v);
            break;
        }
        if(bt != VariantBufferMode::buffer_all)
        {
            if (chr != "" && v.chr != chr)
            {
                // put last variant back
                source.enqueue(v);
                break;
            }
            chr = v.chr;
        }

        if (end >= 0 && ((chr != "" && v.chr != chr) || v.pos > end))
        {
            source.enqueue(v);
            // stop at end -> no more variants
            advance_success = false;
            break;
        }
        // block mode: must not just have no-calls
        if(end == -2 && any_calls && ((chr != "" && v.chr != chr) || (last_end >= 0 && v.pos >= last_end + param)))
        {
            source.enqueue(v);
            break;
        }

        bool all_homref = false;
        if (bt == VariantBufferMode::buffer_block)
        {
            all_homref = true;
            for (Call const & c : v.calls)
            {
                if(!c.isNocall())
                {
                    any_calls = true;
#ifdef DEBUG_VARIANTPROCESSOR
                    std::cerr << "\tfound Call at " << v.chr << ":" << v.pos << " - " << c << "\n";
#endif
                    break;
                }
            }
            for (Call const & c : v.calls)
            {
                if(!(c.isHomref() || c.isNocall()))
                {
                    all_homref = false;
                    break;
                }
            }
            if(all_homref)
            {
                for(auto const & l : v.ambiguous_alleles)
                {
                    if(!l.empty())
                    {
                        all_homref = false;
#ifdef DEBUG_VARIANTPROCESSOR
                        std::cerr << "\tfound ambiguous_alleles at " << v.chr << ":" << v.pos << "\n";
#endif
                        any_calls = true;
                        break;
                    }
                }
            }
        }
        if(!all_homref)
        {
            last_end = std::max(last_end, v.pos + v.len - 1);
        }

        chr = v.chr;

#ifdef DEBUG_VARIANTPROCESSOR
        std::cerr << "\t adding " << v << "\n";
#endif
        add(v);
        ++count;
    }
#ifdef DEBUG_VARIANTPROCESSOR
    std::cerr << "\tend of block (variants still present: " << advance_success << ")" << "\n";
#endif
    return advance_success;
}


/** enqueue a RefVar */
void AbstractVariantProcessingStep::add_variant(int sample, const char * chr, RefVar const & rv, bool het)
{
    Variants vs;
    vs.chr = chr;
    vs.pos = rv.start;
    vs.len = rv.end - rv.start + 1;
    vs.calls.resize(sample + 1);
    vs.variation.push_back(rv);

    vs.calls[sample].ngt = 2;
    if(het)
    {
        vs.calls[sample].gt[0] = 0;
        vs.calls[sample].gt[1] = 1;
    }
    else
    {
        vs.calls[sample].gt[0] = 1;
        vs.calls[sample].gt[1] = 1;
    }
    add(vs);
}

/** enqueue homref block */
void AbstractVariantProcessingStep::add_homref(int sample, const char * chr, int64_t start, int64_t end, bool het)
{
    Variants vs;
    vs.chr = chr;
    vs.pos = start;
    vs.len = end - start + 1;
    vs.calls.resize(sample + 1);
    if (!het)
    {
        vs.calls[sample].ngt = 2;
        vs.calls[sample].gt[0] = 0;
        vs.calls[sample].gt[1] = 0;
    }
    else
    {
        vs.calls[sample].ngt = 1;
        vs.calls[sample].gt[0] = 0;
    }
    add(vs);
}

struct VariantProcessor::VariantProcessorImpl
{
    VariantProcessorImpl() : source(NULL) {}
    VariantProcessorImpl(VariantProcessorImpl const & rhs) :
        mode(rhs.mode),
        param(rhs.param),
        source(rhs.source),
        processing_steps(rhs.processing_steps),
        output_queue(rhs.output_queue) {}

    VariantBufferMode mode;
    int64_t param;
    VariantReader * source;
    std::list< AbstractVariantProcessingStep * > processing_steps;
    std::list<Variants> output_queue;
};

VariantProcessor::VariantProcessor()
{
    _impl = new VariantProcessorImpl();
}

VariantProcessor::VariantProcessor(VariantProcessor const & rhs)
{
    _impl = new VariantProcessorImpl(*rhs._impl);
}

VariantProcessor::~VariantProcessor()
{
    delete _impl;
}

VariantProcessor const & VariantProcessor::operator=(VariantProcessor const & rhs)
{
    if (&rhs == this)
    {
        return *this;
    }
    delete _impl;
    _impl = new VariantProcessorImpl(*rhs._impl);
    return *this;
}


/** set up processing */
void VariantProcessor::addStep(AbstractVariantProcessingStep & step, bool prepend)
{
    if(!prepend)
    {
        _impl->processing_steps.push_back(&step);
    }
    else
    {
        _impl->processing_steps.push_front(&step);
    }
}

/** process a Variant Reader */
void VariantProcessor::setReader(VariantReader & input, VariantBufferMode mode, int64_t param)
{
    _impl->mode = mode;
    _impl->param = param;
    _impl->source = &input;
}

/**
 * @brief Rewind / set region to read
 *
 * @param chr chromosome/contig name
 * @param startpos start position on chr (or -1 for entire chr)
 *
 *
 * Note that this involves clearing internal buffers, which might be slow if
 * the lots of variants are in the pipeline still.
 */
void VariantProcessor::rewind(const char * chr, int64_t startpos)
{
    _impl->output_queue.clear();
    if (!_impl->processing_steps.empty())
    {
        // flush starting at the end
        for (auto & step : boost::adaptors::reverse(_impl->processing_steps))
        {
            step->flush();
        }
    }
    if (_impl->source)
    {
        _impl->source->rewind(chr, startpos);
    }
}

/** Variant output **/
/**
 * @brief Return variant block at current position
 **/
Variants & VariantProcessor::current()
{
    if(!_impl->output_queue.empty())
    {
        return _impl->output_queue.front();
    }
    if (_impl->processing_steps.empty())
    {
        if (!_impl->source)
        {
            error("No input reader and no processing steps.");
        }
#ifdef DEBUG_VARIANTPROCESSOR
        std::cerr << "\t returning " << _impl->source->current() << "\n";
#endif
        return _impl->source->current();
    }
#ifdef DEBUG_VARIANTPROCESSOR
    std::cerr << "\t returning " << _impl->processing_steps.back()->current() << "\n";
#endif
    return _impl->processing_steps.back()->current();
}

/**
 * @brief Advance one line
 * @return true if a variant was retrieved, false otherwise
 */
bool VariantProcessor::advance()
{
    if(!_impl->output_queue.empty())
    {
        _impl->output_queue.pop_front();
        if(!_impl->output_queue.empty())
        {
            return true;
        }
    }
    if (_impl->processing_steps.empty())
    {
        if (!_impl->source)
        {
            error("No input reader and no processing steps.");
        }
        return _impl->source->advance();
    }

    // last step still has things buffered?
    if (_impl->processing_steps.back()->advance())
    {
#ifdef DEBUG_VARIANTPROCESSOR
        std::cerr << "\t advance still has data\n";
#endif
        return true;
    }
    else
    {
#ifdef DEBUG_VARIANTPROCESSOR
        std::cerr << "\t advance buffer refill\n";
#endif
        bool advance_success = false;
        bool any_variants_left = true;
        while(!advance_success && any_variants_left)
        {
            // feed from source if we have one
            if(_impl->source)
            {
                any_variants_left = _impl->processing_steps.front()->add(*(_impl->source), _impl->mode, _impl->param);
            }
            else
            {
                any_variants_left = false;
            }
            auto pstep = _impl->processing_steps.begin();
            auto previous_step = pstep;

#ifdef DEBUG_VARIANTPROCESSOR
            int step = 1;
#endif
#ifdef DEBUG_VARIANTPROCESSOR_STEPS
            std::cerr << "Starting to buffer..." << "\n";
            int istep = 1;
#endif
            while((++pstep) != _impl->processing_steps.end())
            {

#ifdef DEBUG_VARIANTPROCESSOR_STEPS
                std::cerr << "Starting step " << istep << "\n";
#endif
                while((*previous_step)->advance() == true)
                {
#ifdef DEBUG_VARIANTPROCESSOR_STEPS
                    std::cerr << "Adding " << (*previous_step)->current() << "\n";
#endif
                    (*pstep)->add((*previous_step)->current());
#ifdef DEBUG_VARIANTPROCESSOR
                    std::cerr << "\t advancing step " << step << " / success: " << advance_success << "\n";
#endif
                }
#ifdef DEBUG_VARIANTPROCESSOR_STEPS
                std::cerr << "Finished step " << istep++ << "\n";
#endif
#ifdef DEBUG_VARIANTPROCESSOR
                ++step;
#endif
                previous_step = pstep;
            }
            // advance final step
            advance_success = (*previous_step)->advance();
#ifdef DEBUG_VARIANTPROCESSOR
            std::cerr << "\t advance after refill: " << advance_success << "\n";
#endif
        }
        return advance_success;
    }
}

/** put back a variant to the last stage */
void VariantProcessor::putBack(Variants const & v)
{
    _impl->output_queue.push_front(v);
}


/**
 * @brief Remove unused alleles
 */
void trimAlleles(Variants & vars)
{
    // maximally 64 alleles
    uint64_t used = 0;

    if(vars.variation.size() > 64)
    {
        error("trimAlleles only works for up to 64 alt alleles.");
    }

    for(auto & c : vars.calls)
    {
        for(size_t igt = 0; igt < c.ngt; ++igt)
        {
            if (c.gt[igt] > 0)
            {
                used |= (((uint64_t)1) << int8_t(c.gt[igt] - 1));
            }
        }
    }

    for (auto const & l : vars.ambiguous_alleles)
    {
        for(int i : l) {
            if (i > 0)
            {
                used |= (((uint64_t)1) << int8_t(i - 1));
            }
        }
    }

    uint64_t all = (((uint64_t)1) << vars.variation.size()) - 1;
    uint64_t unused = all ^ used;

    if(!unused)
    {
        return;
    }

    int64_t shift = 0;

    size_t al_remap[64];
    size_t nal = vars.variation.size();
    for (size_t pos = 0; pos < nal; ++pos)
    {
        if(unused & 1)
        {
            auto iterator = std::next( vars.variation.begin(), pos + shift );
            vars.variation.erase(iterator);
            shift -= 1;
        }

        al_remap[pos] = pos + shift;
        unused >>= 1;
    }

    for(auto & c : vars.calls)
    {
        for(size_t igt = 0; igt < c.ngt; ++igt)
        {
            if (c.gt[igt] > 0)
            {
                c.gt[igt] = al_remap[c.gt[igt]-1] + 1;
            }
        }
    }
}


}


