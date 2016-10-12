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
 * Table of values
 *
 * \file Table.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef HAPLOTYPES_TABLE_HH
#define HAPLOTYPES_TABLE_HH

#include <memory>
#include <iostream>

namespace table
{
    class Table
    {
    public:
        Table();
        Table(Table const &);
        Table(Table &&);
        ~Table();

        Table &operator=(Table const &);
        Table &operator=(Table &&);

        void set(std::string const & row,
                 std::string const & column,
                 double value
        );

        void set(std::string const & row,
                 std::string const & column,
                 std::string const & value);

        void unset(std::string const & row,
                   std::string const & column);

        void dropRowsWithMissing(std::string const & column);

        bool hasRow(std::string const & row) const;

        std::string getString(std::string const & row,
                              std::string const & column,
                              const char * _def = ".") const;
        double getDouble(std::string const & row,
                         std::string const & column,
                         double _def = std::numeric_limits<double>::quiet_NaN()) const;

    private:
        struct TableImpl;
        std::unique_ptr<TableImpl> _impl;
        friend std::ostream & operator<< (std::ostream & o, Table const & t);
    };

}
#endif //HAPLOTYPES_TABLE_HH
