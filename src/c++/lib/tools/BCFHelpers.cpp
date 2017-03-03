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
 *  \brief Implementation of BCF file format helpers
 *
 * \file BCFHelpers.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "helpers/BCFHelpers.hh"

#include "Error.hh"

#include <cstdio>
#include <sstream>
#include <htslib/vcf.h>
#include <memory>
#include <limits>
#include <set>
#include <boost/algorithm/string.hpp>

/**
 * @brief Helper to get out GT fields
 */

namespace bcfhelpers
{

    namespace _impl
    {
        /** C++-ified version of bcf_get_format_values */
        template <typename target_type_t>
        struct bcf_get_numeric_format
        {
            bcf_get_numeric_format() {}
            /* return true when successful
             */
            void operator()(const bcf_hdr_t *hdr, bcf1_t *line,
                                       const char *tag, int isample,
                                       std::vector<target_type_t> & dest) const
            {
                dest.clear();
                int i;
                int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);

                if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) )
                {
                    return;
                }

                int nsmpl = bcf_hdr_nsamples(hdr);
                if (isample >= nsmpl)
                {
                    return;
                }

                bcf_unpack(line, BCF_UN_FMT);

                for (i = 0; i < line->n_fmt; i++)
                {
                    if ( line->d.fmt[i].id == tag_id )
                    {
                        break;
                    }
                }

                if ( i == line->n_fmt )
                {
                    return;
                }

                bcf_fmt_t *fmt = &line->d.fmt[i];
                int type = fmt->type;
                dest.resize((unsigned long) fmt->n);

                for (i = 0; i < fmt->n; ++i)
                {
                    if ( type == BCF_BT_FLOAT )
                    {
                        static const auto make_missing_float = []() -> float
                        {
                            float f;
                            int32_t _mfloat = 0x7F800001;
                            memcpy(&f, &_mfloat, 4);
                            return f;
                        };
                        static const float bcf_missing_float = make_missing_float();
                        float res = ((float*)(fmt->p + isample*fmt->size))[i];
                        if(res != bcf_missing_float)
                        {
                            dest[i] = target_type_t(res);
                        }
                    }
                    else if (type == BCF_BT_INT8)
                    {
                        int8_t r = ((int8_t*)(fmt->p + isample*fmt->size ))[i];
                        if(r == bcf_int8_vector_end)
                        {
                            break;
                        }
                        if(r != bcf_int8_missing)
                        {
                            dest[i] = target_type_t(r);
                        }
                    }
                    else if (type == BCF_BT_INT16)
                    {
                        int16_t r = (((int16_t*)(fmt->p)) + isample*fmt->size)[i];
                        if(r == bcf_int16_vector_end)
                        {
                            break;
                        }
                        if(r != bcf_int16_missing)
                        {
                            dest[i] = target_type_t(r);
                        }
                    }
                    else if (type == BCF_BT_INT32)
                    {
                        int32_t r = ((int32_t*)(fmt->p + isample*fmt->size))[i];
                        if(r == bcf_int32_vector_end)
                        {
                            break;
                        }
                        if(r != bcf_int32_missing)
                        {
                            dest[i] = target_type_t(r);
                        }
                    }
                    else
                    {
                        std::cerr << "[W] string format field ignored when looking for numeric formats!" << "\n";
                        dest[i] = -1;
                    }
                }
            }
        };

        template <typename type_t>
        struct bcf_get_gts {

            bcf_get_gts() {}

            /**
             * @brief Extract GT from Format field
             *
             */
            void operator() (bcf_fmt_t* gt, int i, int* igt, bool & phased) const
            {
                phased = false;
                type_t *p = (type_t*) (gt->p + i*gt->size);
                int ial;
                for (ial=0; ial < MAX_GT; ial++)
                {
                    igt[ial] = -1;
                }
                for (ial=0; ial < gt->n; ial++)
                {
                    if ( p[ial]==vector_end ) break; /* smaller ploidy */
                    if ( !(p[ial]>>1) || p[ial] == missing ) continue; /* missing allele */
                    int al = (p[ial]>>1)-1;
                    igt[ial] = al;
                    phased = phased || ((p[ial] & 1) != 0);
                }
            }

            static const type_t missing;
            static const type_t vector_end;
        };
        template<> const int8_t bcf_get_gts<int8_t>::missing (bcf_int8_missing);
        template<> const int8_t bcf_get_gts<int8_t>::vector_end (bcf_int8_vector_end);
        template<> const int16_t bcf_get_gts<int16_t>::missing (bcf_int16_missing);
        template<> const int16_t bcf_get_gts<int16_t>::vector_end (bcf_int16_vector_end);
        template<> const int32_t bcf_get_gts<int32_t>::missing (bcf_int32_missing);
        template<> const int32_t bcf_get_gts<int32_t>::vector_end (bcf_int32_vector_end);


        template<typename type_t>
        struct bcf_get_info
        {
            bcf_get_info() {}
            type_t operator()(bcf_info_t * field, int which = 0) const;
        };

        template<>
        int bcf_get_info<int>::operator()(bcf_info_t * field, int which) const
        {
            if(field->type != BCF_BT_CHAR && field->len <= which)
            {
                error("Cannot extract int from non-scalar INFO field (len = %i, requested: %i).",
                      field->len, which);
            }
            switch(field->type)
            {
                case BCF_BT_NULL:
                    return 0;
                case BCF_BT_INT8:
                    return (int) *(((int8_t *)field->vptr) + which );
                case BCF_BT_INT16:
                    return (int) *(((int16_t *)field->vptr) + which );
                case BCF_BT_INT32:
                    return (int) *(((int32_t *)field->vptr) + which );
                case BCF_BT_FLOAT:
                    return (int) *(((float*)field->vptr) + which );
                case BCF_BT_CHAR:
                    if(which > 0)
                    {
                        error("Cannot extract int %i from string INFO field", which);
                    }
                    return atoi((const char *)field->vptr);
                default:
                    break;
            }
            return -1;
        }

        template<>
        float bcf_get_info<float>::operator()(bcf_info_t * field, int which) const
        {
            switch(field->type)
            {
                case BCF_BT_NULL:
                    return std::numeric_limits<float>::quiet_NaN();
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                {
                    const bcf_get_info<int> gi;
                    return (float)gi(field, which);
                }
                case BCF_BT_FLOAT:
                    if(field->len <= which)
                    {
                        error("Cannot extract int from non-scalar INFO field (len = %i, requested: %i).",
                              field->len, which);

                    }
                    return *(((float*)field->vptr) + which );
                case BCF_BT_CHAR:
                    if(which > 0)
                    {
                        error("Cannot extract int %i from string INFO field", which);
                    }
                    return (float) atof((const char *)field->vptr);
                default:break;
            }
            return std::numeric_limits<float>::quiet_NaN();
        }

        template<>
        std::string bcf_get_info<std::string>::operator()(bcf_info_t * field, int) const
        {
            char num[256];
            switch(field->type)
            {
                case BCF_BT_NULL:
                    return std::string("NULL");
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                    snprintf(num, 256, "%i", field->v1.i);
                    return std::string(num);
                case BCF_BT_FLOAT:
                    snprintf(num, 256, "%g", field->v1.f);
                    return std::string(num);
                case BCF_BT_CHAR:
                    if(field->vptr && field->len > 0)
                    {
                        return std::string((const char *)field->vptr, (unsigned long) field->len);
                    }
                    break;
                default:break;
            }
            return std::string();
        }

    }

    void bcfHeaderHG19(bcf_hdr_t * header)
    {
        bcf_hdr_append(header, "##reference=hg19");
        bcf_hdr_append(header, "##contig=<ID=chr1,length=249250621>");
        bcf_hdr_append(header, "##contig=<ID=chr2,length=243199373>");
        bcf_hdr_append(header, "##contig=<ID=chr3,length=198022430>");
        bcf_hdr_append(header, "##contig=<ID=chr4,length=191154276>");
        bcf_hdr_append(header, "##contig=<ID=chr5,length=180915260>");
        bcf_hdr_append(header, "##contig=<ID=chr6,length=171115067>");
        bcf_hdr_append(header, "##contig=<ID=chr7,length=159138663>");
        bcf_hdr_append(header, "##contig=<ID=chr8,length=146364022>");
        bcf_hdr_append(header, "##contig=<ID=chr9,length=141213431>");
        bcf_hdr_append(header, "##contig=<ID=chr10,length=135534747>");
        bcf_hdr_append(header, "##contig=<ID=chr11,length=135006516>");
        bcf_hdr_append(header, "##contig=<ID=chr12,length=133851895>");
        bcf_hdr_append(header, "##contig=<ID=chr13,length=115169878>");
        bcf_hdr_append(header, "##contig=<ID=chr14,length=107349540>");
        bcf_hdr_append(header, "##contig=<ID=chr15,length=102531392>");
        bcf_hdr_append(header, "##contig=<ID=chr16,length=90354753>");
        bcf_hdr_append(header, "##contig=<ID=chr17,length=81195210>");
        bcf_hdr_append(header, "##contig=<ID=chr18,length=78077248>");
        bcf_hdr_append(header, "##contig=<ID=chr19,length=59128983>");
        bcf_hdr_append(header, "##contig=<ID=chr20,length=63025520>");
        bcf_hdr_append(header, "##contig=<ID=chr21,length=48129895>");
        bcf_hdr_append(header, "##contig=<ID=chr22,length=51304566>");
        bcf_hdr_append(header, "##contig=<ID=chrX,length=155270560>");
        bcf_hdr_append(header, "##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">");
        bcf_hdr_append(header, "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that could not be imported.\">");
        bcf_hdr_append(header, "##FORMAT=<ID=AGT,Number=1,Type=String,Description=\"Genotypes at ambiguous locations\">");
        bcf_hdr_append(header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(header, "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
        bcf_hdr_append(header, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
        bcf_hdr_append(header, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">");
        bcf_hdr_append(header, "##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">");
    }


    /** extract chrom / pos / length */
    void getLocation(bcf_hdr_t * hdr, bcf1_t * rec, int64_t & refstart, int64_t & refend)
    {
        bcf_unpack(rec, BCF_UN_STR);
        refstart = rec->pos;
        refend = refstart;

        int endfield = getInfoInt(hdr, rec, "END", -1);

        if(endfield > 0)
        {
            // if there is an end field, don't validate the ref allele
            refend = endfield-1;
        }
        else if(rec->d.allele)
        {
            refend = refstart + strlen(rec->d.allele[0]) - 1;
            if (strchr(rec->d.allele[0], '.') || strchr(rec->d.allele[0], '-'))
            {
                // length might be inaccurate now
                refend = refstart;
                throw bcfhelpers::importexception(
                    std::string("[W] Unsupported REF allele with undefined length: ") +
                    std::string(rec->d.allele[0]));
            }
        }
    }

    /**
     * @brief Retrieve an info field as an integer
     *
     * @param result the default to return if the field is not present
     */
    std::string getInfoString(bcf_hdr_t * header, bcf1_t * line, const char * field, const char * def_result)
    {
        std::string result = def_result;
        bcf_info_t * info_ptr  = bcf_get_info(header, line, field);
        if(info_ptr)
        {
            static const bcfhelpers::_impl::bcf_get_info<std::string> i;
            result = i(info_ptr);
        }
        return result;
    }

    /**
     * @brief Retrieve an info field as an integer
     *
     * @param result the default to return if the field is not present
     */
    int getInfoInt(bcf_hdr_t * header, bcf1_t * line, const char * field, int result)
    {
        bcf_unpack(line, BCF_UN_INFO);
        bcf_info_t * info_ptr  = bcf_get_info(header, line, field);
        if(info_ptr)
        {
            static const bcfhelpers::_impl::bcf_get_info<int> i;
            result = i(info_ptr);
        }
        return result;
    }

    std::vector<int> getInfoInts(bcf_hdr_t * header, bcf1_t * line, const char * field)
    {
        std::vector<int> result;
        bcf_unpack(line, BCF_UN_INFO);
        bcf_info_t * info_ptr  = bcf_get_info(header, line, field);
        if(info_ptr)
        {
            static const bcfhelpers::_impl::bcf_get_info<int> i;
            for(int j = 0; j < info_ptr->len; ++j)
            {
                result.push_back(i(info_ptr, j));
            }
        }
        return result;
    }

    /**
     * @brief Retrieve an info field as a double
     *
     * @return the value or NaN
     */
    float getInfoFloat(bcf_hdr_t * header, bcf1_t * line, const char * field)
    {
        float result = std::numeric_limits<float>::quiet_NaN();
        bcf_unpack(line, BCF_UN_INFO);
        bcf_info_t * info_ptr  = bcf_get_info(header, line, field);
        if(info_ptr)
        {
            static const bcfhelpers::_impl::bcf_get_info<float> i;
            result = i(info_ptr);
        }
        return result;
    }

    std::vector<float> getInfoFloats(bcf_hdr_t * header, bcf1_t * line, const char * field)
    {
        std::vector<float> result;
        bcf_unpack(line, BCF_UN_INFO);
        bcf_info_t * info_ptr  = bcf_get_info(header, line, field);
        if(info_ptr)
        {
            static const bcfhelpers::_impl::bcf_get_info<float> i;
            for(int j = 0; j < info_ptr->len; ++j)
            {
                result.push_back(i(info_ptr, j));
            }
        }
        return result;
    }

    /**
     * @brief Retrieve an info flag
     *
     * @return true of the flag is set
     */
    bool getInfoFlag(bcf_hdr_t * hdr, bcf1_t * line, const char * field)
    {
        bcf_unpack(line, BCF_UN_INFO);
        return bcf_get_info_flag(hdr, line, field, nullptr, 0) == 1;
    }

    /**
     * @brief Read the GT field
     */
    void getGT(bcf_hdr_t * header, bcf1_t * line, int isample, int * gt, int & ngt, bool & phased)
    {
        bcf_unpack(line, BCF_UN_FMT);
        bcf_fmt_t * fmt_ptr = bcf_get_fmt(header, line, "GT");
        if (fmt_ptr)
        {
            ngt = fmt_ptr->n;
            if(ngt > MAX_GT)
            {
                std::cerr << "[W] Found a variant with more " << ngt << " > " << MAX_GT << " (max) alt alleles. These become no-calls." << "\n";
                for (int i = 0; i < MAX_GT; ++i)
                {
                    gt[i] = -1;
                }
                phased = false;
                return;
            }

            // check if all alleles are populated
            switch (fmt_ptr->type) {
                case BCF_BT_INT8:  {static const bcfhelpers::_impl::bcf_get_gts<int8_t> b; b(fmt_ptr, isample, gt, phased);} break;
                case BCF_BT_INT16: {static const bcfhelpers::_impl::bcf_get_gts<int16_t> b; b(fmt_ptr, isample, gt, phased);} break;
                case BCF_BT_INT32: {static const bcfhelpers::_impl::bcf_get_gts<int32_t> b; b(fmt_ptr, isample, gt, phased);} break;
                default: error("Unsupported GT type: %d at %s:%d\n", fmt_ptr->type,
                               header->id[BCF_DT_CTG][line->rid].key, line->pos+1);  break;
            }
        }
        else
        {
            ngt = 0;
            phased = false;
        }
    }

    /** read GQ(X) -- will use in this order: GQX, GQ, -1
     *
     * This is somewhat Illumina-specific, TODO: refactor to avoid using this function in a general setting
     * */
    void getGQ(const bcf_hdr_t * header, bcf1_t * line, int isample, float & gq)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<float> gf;

        std::vector<float> values;
        gf(header, line, "GQX", isample, values);
        if(values.empty())
        {
            gf(header, line, "GQ", isample, values);
        }
        if(values.size() > 1)
        {
            std::cerr << "[W] too many GQ fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos << "\n";
        }
        if(values.empty())
        {
            gq = -1;
        }
        else
        {
            gq = values[0];
        }
    }

    /** read AD */
    void getAD(const bcf_hdr_t * header, bcf1_t * line, int isample, int *ad, int max_ad)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;

        std::vector<int> values;
        gf(header, line, "AD", isample, values);
        if(max_ad < (int)values.size())
        {
            std::cerr << "[W] too many AD fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
                      << " max_ad = " << max_ad << " retrieved: " << values.size()
                      << "\n";
        }
        for(size_t q = 0; q < values.size(); ++q)
        {
            ad[q] = values[q];
        }
    }

    /** read DP(I) -- will use in this order: DP, DPI, -1 */
    void getDP(const bcf_hdr_t * header, bcf1_t * line, int isample, int & dp)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;

        std::vector<int> values;
        gf(header, line, "DP", isample, values);
        if(values.empty())
        {
            gf(header, line, "DPI", isample, values);
        }
        if(values.size() > 1)
        {
            std::cerr << "[W] too many DP fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos << "\n";
        }
        if(values.empty())
        {
            dp = 0;
        }
        else
        {
            dp = values[0];
        }
    }

    /** read a format field as a single int. */
    int getFormatInt(const bcf_hdr_t * header, bcf1_t * line, const char * field, int isample, int defaultresult)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;

        std::vector<int> values;
        gf(header, line, field, isample, values);
        if(values.size() > 1)
        {
            std::ostringstream os;
            os << "[W] too many " << field << " fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos;
            throw importexception(os.str());
        }
        if(values.size() == 1)
        {
            return values[0];
        }
        return defaultresult;
    }

    std::vector<int> getFormatInts(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;
        std::vector<int> result;
        gf(header, line, field, isample, result);
        return result;
    }

    /** read a format field as a single float. default return value is NaN */
    float getFormatFloat(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
    {
        using namespace _impl;
        float result = std::numeric_limits<float>::quiet_NaN();
        static const bcf_get_numeric_format<float> gf;

        std::vector<float> values;
        gf(header, line, field, isample, values);
        if(values.size() > 1)
        {
            std::ostringstream os;
            os << "[W] too many " << field << " fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos;
            throw importexception(os.str());
        }
        if(values.size() == 1)
        {
            result = values[0];
        }
        return result;
    }

    std::vector<float> getFormatFloats(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<float> gf;
        std::vector<float> result;
        gf(header, line, field, isample, result);
        return result;
    }

    /** read a format field as a single double. result will not be overwritten on failure */
    std::string getFormatString(const bcf_hdr_t * hdr, bcf1_t * line, const char * field, int isample, const char * result)
    {
        int nsmpl = bcf_hdr_nsamples(hdr);
        if (isample >= nsmpl)
        {
            return result;
        }

        bcf_unpack(line, BCF_UN_FMT);
        int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, field);

        if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) )
        {
            return result;
        }

        if ( !(line->unpacked & BCF_UN_FMT) )
        {
            bcf_unpack(line, BCF_UN_FMT);
        }

        // index in format fields
        int i = 0;
        for (i = 0; i < line->n_fmt; i++)
        {
            if ( line->d.fmt[i].id == tag_id )
            {
                break;
            }
        }

        if ( i == line->n_fmt )
        {
            return result;
        }

        bcf_fmt_t *fmt = &line->d.fmt[i];

        if(fmt == NULL)
        {
            return result;
        }

        int type = fmt->type;

        if(fmt->n < 1)
        {
            return result;
        }

        std::string str_result = result;
        if ( type == BCF_BT_FLOAT )
        {
            static const auto make_missing_float = []() -> float
            {
                float f;
                int32_t _mfloat = 0x7F800001;
                memcpy(&f, &_mfloat, 4);
                return f;
            };
            static const float bcf_missing_float = make_missing_float();
            float res = *((float*)(fmt->p + isample*fmt->size));
            if(res != bcf_missing_float)
            {
                str_result = std::to_string(res);
            }
        }
        else if (type == BCF_BT_INT8)
        {
            int8_t r = *((int8_t*)(fmt->p + isample*fmt->size ));
            if(r != bcf_int8_missing && r != bcf_int8_vector_end)
            {
                str_result = std::to_string(r);
            }
        }
        else if (type == BCF_BT_INT16)
        {
            int16_t r = *(((int16_t*)(fmt->p)) + isample*fmt->size);
            if(r != bcf_int16_missing && r != bcf_int16_vector_end)
            {
                str_result = std::to_string(r);
            }
        }
        else if (type == BCF_BT_INT32)
        {
            int32_t r = *((int32_t*)(fmt->p + isample*fmt->size));
            if(r != bcf_int32_missing && r == bcf_int32_vector_end)
            {
                str_result = std::to_string(r);
            }
        }
        else
        {
            const char * src = (const char *)fmt->p + isample*fmt->size;
            if(src && fmt->size)
            {
                str_result = std::string(src, (unsigned long) fmt->size);
                // deal with 0 padding
                str_result.resize(strlen(str_result.c_str()));
            }
        }

        return str_result;
    }

    /** update format string for a single sample.  */
    void setFormatStrings(const bcf_hdr_t * hdr, bcf1_t * line, const char * field,
                         const std::vector<std::string> & formats)
    {
        // TODO this can probably be done faster / better
        std::unique_ptr<const char *[]> p_fmts = std::unique_ptr<const char *[]>(new const char *[line->n_sample]);
        bcf_unpack(line, BCF_UN_FMT);
        bool any_nonempty = false;
        if(formats.size() != line->n_sample)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field
               << " " << hdr->id[BCF_DT_ID][line->rid].key << ":" << line->pos
               << " -- we have " << formats.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
        for (int si = 0; si < line->n_sample; ++si)
        {
            p_fmts.get()[si] = formats[si].c_str();
            if(!formats[si].empty())
            {
                any_nonempty = true;
            }
        }
        int res;
        if(any_nonempty)
        {
            res = bcf_update_format_string(hdr, line, field, p_fmts.get(), line->n_sample);
        }
        else
        {
            res = bcf_update_format_string(hdr, line, field, NULL, 0);
        }

        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << hdr->id[BCF_DT_ID][line->rid].key << ":" << line->pos
               << " -- we have " << formats.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
    }

    /** update format with single float values.  */
    void setFormatFloats(const bcf_hdr_t * header, bcf1_t * line, const char * field,
                         const std::vector<float> & value)
    {
        if(value.empty())
        {
            int res = bcf_update_format(header, line, field, NULL, 0, 0);
            if(res != 0)
            {
                std::ostringstream os;
                os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
                   << " -- we have " << value.size() << " for " << line->n_sample << " samples";

                throw importexception(os.str());
            }
            return;
        }
        std::unique_ptr<float[]> p_dbl = std::unique_ptr<float[]>(new float[line->n_sample]);

        for(size_t i = 0; i < line->n_sample; ++i)
        {
            if(i < value.size())
            {
                p_dbl.get()[i] = value[i];
            }
            else
            {
                p_dbl.get()[i] = std::numeric_limits<float>::quiet_NaN();
            }
        }

        int res = bcf_update_format_float(header, line, field, p_dbl.get(), line->n_sample);
        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
               << " -- we have " << value.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
    }

    void setFormatInts(const bcf_hdr_t * header, bcf1_t * line, const char * field,
                       const std::vector<int> & value)
    {
        if(value.empty())
        {
            int res = bcf_update_format(header, line, field, NULL, 0, 0);
            if(res != 0)
            {
                std::ostringstream os;
                os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
                   << " -- we have " << value.size() << " for " << line->n_sample << " samples";

                throw importexception(os.str());
            }
            return;
        }
        std::unique_ptr<int[]> p_dbl = std::unique_ptr<int[]>(new int[line->n_sample]);

        for(size_t i = 0; i < line->n_sample; ++i)
        {
            if(i < value.size())
            {
                p_dbl.get()[i] = value[i];
            }
            else
            {
                p_dbl.get()[i] = bcf_int32_missing;
            }
        }

        int res = bcf_update_format_int32(header, line, field, p_dbl.get(), line->n_sample);
        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
               << " -- we have " << value.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
    }

    /** return sample names from header */
    std::list<std::string> getSampleNames(const bcf_hdr_t * hdr)
    {
        std::list<std::string> l;
        for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
        {
            std::string samplename = hdr->samples[i];
            if(samplename == "*")
            {
                std::cerr << "Skipping sample named '*'" << "\n";
                continue;
            }
            l.push_back(samplename);
        }
        return l;
    }

    /** return number of reference padding bases */
    int isRefPadded(bcf1_t * line)
    {
        bcf_unpack(line, BCF_UN_SHR);
        if(line->n_allele == 1)
        {
            return 0;
        }

        if(!line->d.allele)
        {
            throw bcfhelpers::importexception("Empty REF allele -- this is probably because parsing this record failed.");
        }

        const char * ref = line->d.allele[0];
        const int reflen = (int) strlen(ref);

        int max_match = reflen;
        for(int al = 1; al < line->n_allele; ++al)
        {
            const char * alt = line->d.allele[al];
            // symbolic or missing ALT
            if(strcmp(alt, ".") == 0 || *alt == '<')
            {
                max_match = 0;
                break;
            }
            int rpos = 0;
            for(rpos = 0; rpos < reflen; ++rpos, ++alt)
            {
                const char rb = *(ref + rpos);
                if(*alt == 0 || *alt != rb)
                {
                    break;
                }
            }
            max_match = std::min(rpos, max_match);
        }
        return max_match;
    }

    std::pair<AlleleType, std::string> classifyAlleleString(std::string allele_string)
    {
        boost::to_upper(allele_string);

        if (allele_string == "<DEL>" || allele_string == "*")
        {
            return std::make_pair(AlleleType::SYMBOLIC_DEL, std::string());
        }
        else if (allele_string == "" || allele_string == ".")
        {
            return std::make_pair(AlleleType::MISSING, std::string("."));
        }
        else if(allele_string.find('<') != std::string::npos || allele_string.find('>') != std::string::npos)
        {
            return std::make_pair(AlleleType::SYMBOLIC_OTHER, allele_string);
        }
        else
        {
            // check alleles
            for (size_t qq = 0; qq < allele_string.size(); ++qq)
            {
                /* A	A	Adenine */
                /* C	C	Cytosine */
                /* G	G	Guanine */
                /* T	T	Thymine */
                /* U	U	Uracil */
                /* R	A or G	puRine */
                /* Y	C, T or U	pYrimidines */
                /* K	G, T or U	bases which are Ketones */
                /* M	A or C	bases with aMino groups */
                /* S	C or G	Strong interaction */
                /* W	A, T or U	Weak interaction */
                /* B	not A (i.e. C, G, T or U)	B comes after A */
                /* D	not C (i.e. A, G, T or U)	D comes after C */
                /* H	not G (i.e., A, C, T or U)	H comes after G */
                /* V	neither T nor U (i.e. A, C or G)	V comes after U */
                /* N	A C G T U	Nucleic acid */
                /* X	masked */
                /* -/.	gap of indeterminate length - not supported */
                if (allele_string[qq] != 'A' && \
                   allele_string[qq] != 'C' && \
                   allele_string[qq] != 'G' && \
                   allele_string[qq] != 'T' && \
                   allele_string[qq] != 'U' && \
                   allele_string[qq] != 'R' && \
                   allele_string[qq] != 'Y' && \
                   allele_string[qq] != 'K' && \
                   allele_string[qq] != 'M' && \
                   allele_string[qq] != 'S' && \
                   allele_string[qq] != 'W' && \
                   allele_string[qq] != 'B' && \
                   allele_string[qq] != 'D' && \
                   allele_string[qq] != 'H' && \
                   allele_string[qq] != 'V' && \
                   allele_string[qq] != 'N' && \
                   allele_string[qq] != 'X')
                {
                    return std::make_pair(AlleleType::UNKNOWN, allele_string);
                }
            }
            return std::make_pair(AlleleType::NUC, allele_string);
        }
    }

    std::ostream & operator<<(std::ostream & o, const AlleleType at)
    {
        switch(at)
        {
            case AlleleType::NUC:
                o << "NUC";
                break;
            case AlleleType::MISSING:
                o << "MISSING";
                break;
            case AlleleType::SYMBOLIC_DEL:
                o << "SYMBOLIC_DEL";
                break;
            case AlleleType::SYMBOLIC_OTHER:
                o << "SYMBOLIC_OTHER";
                break;
            case AlleleType::UNKNOWN:
                o << "UNKNOWN";
                break;
        }
        return o;
    }
} // namespace bcfhelpers
