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
 * 
 * Three-sequence multiple alignment (calling muscle)
 *
 * \file MultipleAlignment.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "MultipleAlignment.hh"

#include <memory>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>

#include <boost/filesystem.hpp>
#ifndef __MACH__
#include <unistd.h>
#else
#include <mach-o/dyld.h>
#endif

#include "Error.hh"
#include "Alignment.hh"
#include "helpers/StringUtil.hh"

// #define DEBUG_MULTIALIGN

/**
 * Find path to current executable
 */
#ifndef __MACH__
std::string executablePath()
{
    char pathbuf[4097];
    ssize_t len;

    if ((len = readlink("/proc/self/exe", pathbuf, 4096)) == -1)
    {
        error("Could not determine my own executable's path.");
    }
    pathbuf[len] = '\0';
    return std::string(pathbuf);
}
#else
std::string executablePath()
{
    char pathbuf[4096];
    uint32_t len = sizeof(pathbuf);
    if (_NSGetExecutablePath(pathbuf, &len) != 0)
    {
        error("Could not determine my own executable's path.");
    }
    return std::string(pathbuf);
}
#endif

struct MultipleAlignmentImpl
{
    MultipleAlignmentImpl()
    {
        boost::filesystem::path p(executablePath());
        boost::filesystem::path tp = p.parent_path() / boost::filesystem::path("muscle");
#ifdef HAS_MUSCLE
        if(!boost::filesystem::exists(tp))
        {
            error("Cannot find bundled version of muscle. I'd expect it here: %s", tp.c_str());
        }        
        path_to_muscle = tp.native();
#else
        path_to_muscle = "";
#endif
    }
    std::string path_to_muscle;
};

MultipleAlignment::MultipleAlignment() 
{
    _impl = new MultipleAlignmentImpl();
}

MultipleAlignment::~MultipleAlignment()
{
    delete _impl;
}


/** run command and capture output */
std::string exec(std::string const & cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe)
    {
        error("Error running command '%s'", cmd.c_str());
    }
    char buffer[1024];
    std::string result = "";
    while(!feof(pipe)) 
    {
        if(fgets(buffer, 1024, pipe) != NULL)
        {
            result += buffer;
        }
    }
    int estatus = pclose(pipe);
    if(WEXITSTATUS(estatus) != 0)
    {
        error("Command '%s' exited with status %i / output %s", cmd.c_str(), estatus, result.c_str());
    }
    return result;
}

/**
 * Compute three-sequence multiple alignment, return strings padded with '-'
 */
void MultipleAlignment::align(
    std::string const& ref,
    std::string const& alt1,
    std::string const& alt2,
    std::string & padded_ref,
    std::string & padded_alt1,
    std::string & padded_alt2)
{
#ifdef NO_MUSCLE
    error("Cannot perform multiple alignment when Muscle3 is not available.");
#endif
#ifndef HAS_MUSCLE
    error("Cannot perform multiple alignment when Muscle3 is not available.");
#endif
    boost::filesystem::path tempin = boost::filesystem::unique_path(boost::filesystem::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.fa");
    boost::filesystem::path tempout = boost::filesystem::unique_path(boost::filesystem::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.fa");
    
    {
        std::ofstream minput(tempin.native());
        minput << ">r" << "\n" << ref << "\n";
        minput << ">a1" << "\n" << alt1 << "\n";
        minput << ">a2" << "\n" << alt2 << "\n";
    }
    
    exec(_impl->path_to_muscle + " -in " + tempin.native() + " -out " + tempout.native() + " -diags -quiet");    
    
    {
        std::ifstream moutput(tempout.native());
        std::string line;
        std::string * current = NULL;
        while (std::getline(moutput, line))
        {
            stringutil::replaceAllInplace(line, "\n", "");
            stringutil::replaceAllInplace(line, "\r", "");
            stringutil::replaceAllInplace(line, " ", "");
            if(line == ">r")
            {
                current = &padded_ref;
            }
            else if(line == ">a1")
            {
                current = &padded_alt1;
            }
            else if(line == ">a2")
            {
                current = &padded_alt2;
            }
            else if(current)
            {
                *current += line;
            }
            else
            {
                error("Could not parse muscle output: %s", line.c_str());
            }
        }
    }

#ifdef DEBUG_MULTIALIGN
    std::cerr << padded_ref << "\n";
    std::cerr << padded_alt1 << "\n";
    std::cerr << padded_alt2 << "\n";
#endif

    if(padded_ref.size() != padded_alt1.size() || padded_ref.size() != padded_alt2.size())
    {
        error("Read multi-aligned sequences of different lengths: %s/%s/%s", padded_ref.c_str(), padded_alt1.c_str(), padded_alt2.c_str());
    }
}
