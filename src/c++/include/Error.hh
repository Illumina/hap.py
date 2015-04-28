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
 *  \brief Error handling helper
 *  
 * \file Error.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cstdio>
#include <cstdarg>
#include <stdexcept>

#ifndef ERRMSG_MAX
#define ERRMSG_MAX 2048
#endif

#ifdef _DEBUG

#ifdef __GNUC__

#include <cstdio>
#include <execinfo.h>
#include <signal.h>
#include <cstdlib>

#endif

inline void error_stop(const char * f, int l, const char *format, ...) {
	fprintf(stderr, "\n[DEBUG] Error at %s:%i\n", f, l);
	{ 	// extra Debug print when errors occur
		va_list ap;
		va_start(ap, format);
		vfprintf(stderr, format, ap);
		va_end(ap);
		fprintf(stderr, "\n");
	}

	// in GNUC, do backtrace
#ifdef __GNUC__
	void *array[20];

	size_t size;
	size = backtrace(array, 20);
	backtrace_symbols_fd(array, size, 2);
#endif
	fprintf(stderr, "\n[DEBUG END]\n");

	char errmsg[ERRMSG_MAX];
	va_list ap;
	va_start(ap, format);
	vsnprintf(errmsg, ERRMSG_MAX, format, ap);
	va_end(ap);
	throw std::runtime_error(errmsg);
}

#ifdef error
#undef error
#endif

#define error(...) do { 	\
	error_stop(__FILE__, __LINE__, __VA_ARGS__); \
} while (0)

#else  // _DEBUG

static inline void error(const char *format, ...)
{
	char errmsg[ERRMSG_MAX];
	va_list ap;
	va_start(ap, format);
	vsnprintf(errmsg, ERRMSG_MAX, format, ap);
	va_end(ap);
	throw std::runtime_error(errmsg);
}

#endif
