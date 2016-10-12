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
 * \file Table.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>
#include <list>
#include <set>

#include "helpers/Table.hh"
#include "Error.hh"

 namespace table
 {
     namespace _tableImpl
     {
         struct OptionalValue
         {
             OptionalValue() :
                 dvalue(std::numeric_limits<double>::quiet_NaN()), status(MISSING) {}

             explicit OptionalValue(double _value) :
                 dvalue(_value), status(DOUBLE) {}

             explicit OptionalValue(std::string const & s) :
                 svalue(s), status(STRING) {}

             OptionalValue & operator=(double val)
             {
                 dvalue = val;
                 svalue = "";
                 status = DOUBLE;
                 return *this;
             }

             OptionalValue & operator=(std::string const & s)
             {
                 dvalue = std::numeric_limits<double>::quiet_NaN();
                 svalue = s;
                 status = STRING;
                 return *this;
             }

             void reset()
             {
                 svalue = "";
                 dvalue = std::numeric_limits<double>::quiet_NaN();
                 status = MISSING;
             }

             enum Status {
                 MISSING,
                 DOUBLE,
                 STRING
             };

             double dvalue;
             std::string svalue;
             Status status;

             std::string toString() const
             {
                 switch(status)
                 {
                     case STRING: return svalue;
                     case DOUBLE: return std::to_string(dvalue);
                     case MISSING: return ".";
                 }
                 return ".";
             }

             operator bool() { return status != MISSING; }
             operator double() {
                 return status == DOUBLE ? dvalue :
                        std::numeric_limits<double>::quiet_NaN();
             }
             operator std::string() {
                 return toString();
             }
         };

         typedef std::unordered_map<std::string, OptionalValue> Column;
     }

     struct Table::TableImpl
     {
         std::unordered_map<std::string, _tableImpl::Column> rows;
     };

     Table::Table() : _impl(new TableImpl()) {}
     Table::Table(Table const & rhs) : _impl(new TableImpl(*rhs._impl)) {}
     Table::Table(Table && rhs) : _impl(std::move(rhs._impl)) {}
     Table::~Table() {}

     Table & Table::operator=(Table const & rhs)
     {
         if(&rhs == this)
         {
             return *this;
         }
         *_impl = *rhs._impl;
         return *this;
     }
     Table & Table::operator=(Table && rhs)
     {
         if(&rhs == this)
         {
             return *this;
         }
         _impl = std::move(rhs._impl);
         return *this;
     }

     void Table::set(std::string const & row,
                     std::string const & column,
                     double value
     )
     {
         auto row_it = _impl->rows.find(row);
         if(row_it == _impl->rows.end())
         {
             row_it = _impl->rows.emplace(row, _tableImpl::Column()).first;
         }
         row_it->second[column] = value;
     }

     void Table::set(std::string const & row,
                     std::string const & column,
                     std::string const & value)
     {
         auto row_it = _impl->rows.find(row);
         if(row_it == _impl->rows.end())
         {
             row_it = _impl->rows.emplace(row, _tableImpl::Column()).first;
         }
         row_it->second[column] = value;
     }

     void Table::unset(std::string const & row,
                std::string const & column)
     {
         auto row_it = _impl->rows.find(row);
         if(row_it == _impl->rows.end())
         {
             return;
         }
         auto col_it = row_it->second.find(column);
         if(col_it == row_it->second.end())
         {
             return;
         }
         row_it->second.erase(col_it);
     }

     bool Table::hasRow(std::string const & row) const
     {
         return _impl->rows.find(row) != _impl->rows.end();
     }

     std::string Table::getString(std::string const & row,
                           std::string const & column,
                           const char * _def) const
     {
         auto row_it = _impl->rows.find(row);
         if(row_it == _impl->rows.end())
         {
             return _def;
         }
         auto col_it = row_it->second.find(column);
         if(col_it == row_it->second.end())
         {
             return _def;
         }
         return col_it->second;
     }

     double Table::getDouble(std::string const & row,
                           std::string const & column,
                      double _def) const
     {
         auto row_it = _impl->rows.find(row);
         if(row_it == _impl->rows.end())
         {
             return _def;
         }
         auto col_it = row_it->second.find(column);
         if(col_it == row_it->second.end())
         {
             return _def;
         }
         return col_it->second;
     }

     void Table::dropRowsWithMissing(std::string const & column)
     {
         typedef std::unordered_map<std::string, _tableImpl::Column>::iterator row_iterator;
         std::list<row_iterator> to_delete;
         for(auto row_it = _impl->rows.begin(); row_it != _impl->rows.end(); ++row_it)
         {
             auto col_it = row_it->second.find(column);
             if(col_it == row_it->second.end())
             {
                 to_delete.push_back(row_it);
             }
             else if(col_it->second.status == _tableImpl::OptionalValue::MISSING)
             {
                 to_delete.push_back(row_it);
             }
         }

         for(auto const & row_it : to_delete)
         {
             _impl->rows.erase(row_it);
         }
     }

     std::ostream & operator<< (std::ostream & o, Table const & t)
     {
         std::set<std::string> columns;
         for(auto const & r : t._impl->rows)
         {
             for(auto const & c : r.second)
             {
                 columns.insert(c.first);
             }
         }

         bool first = true;
         for(auto const & c : columns)
         {
             if(!first)
             {
                 o << "\t";
             }
             else
             {
                 first = false;
             }
             o << c;
         }
         o << "\n";

         for(auto const & r : t._impl->rows)
         {
             first = true;
             for(auto const & c : columns)
             {
                 if(!first)
                 {
                     o << "\t";
                 }
                 else
                 {
                     first = false;
                 }
                 auto c_it = r.second.find(c);
                 if(c_it == r.second.end())
                 {
                     o << ".";
                 }
                 else
                 {
                     o << c_it->second.toString();
                 }
             }
             o << "\n";
         }

         return o;
     }
 }
