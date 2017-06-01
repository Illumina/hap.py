#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt

macro(get_compiler_version compiler_version)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE ${compiler_version})
    STRING(REGEX REPLACE "(\r?\n)+$" "" ${compiler_version} "${${compiler_version}}")
endmacro()

# clang doesn't make finding the version easy for us...
macro(get_clang_version compiler_version)
    execute_process(COMMAND bash -c "${CMAKE_CXX_COMPILER} -v 2>&1 | awk '{printf $3; exit}'" OUTPUT_VARIABLE ${compiler_version})
endmacro()

macro(test_min_compiler compiler_version min_compiler_version compiler_label)
    if (${compiler_version} VERSION_LESS ${min_compiler_version})
        message (FATAL_ERROR "Unsupported ${compiler_label} version: ${compiler_version}: "
                             "only versions >= ${min_compiler_version} are supported")
    endif ()
endmacro()


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_compiler_version(compiler_version)
    test_min_compiler(${compiler_version} "4.8.0" "g++")
    message (STATUS "using compiler: g++ version ${compiler_version}")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    get_clang_version(compiler_version)
    message (STATUS "using compiler: clang++ version ${compiler_version}")
else ()
    message (STATUS "using compiler: ${CMAKE_CXX_COMPILER_ID}")
endif ()


#
# set compile flags, and modify by compiler/version:
#
set (GNU_COMPAT_COMPILER ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") ))

# start with warning flags:
if (GNU_COMPAT_COMPILER)
    set (CXX_WARN_FLAGS "-Wall -Wextra -Wunused -Wpointer-arith -Winit-self -Wredundant-decls -pedantic -Wunused-parameter -Wdisabled-optimization -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wuninitialized")

    if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED")
    endif ()
endif ()

#
# add extra compiler specific flags:
#
if     (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    if (${compiler_version} VERSION_LESS "4.7")
        message(FATAL_ERROR "You need at least GCC 4.7 to compile hap.py")
    endif ()

    set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wlogical-op")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    # add LTO for release versions
    if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    endif ()
    set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wno-unused-function -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wno-unknown-pragmas -Wmissing-prototypes -Wunused-exception-parameter -Wbool-conversion -Wempty-body -Wimplicit-fallthrough -Wsizeof-array-argument -Wstring-conversion -Wno-unused-parameter -Wno-missing-prototypes -Wno-unknown-warning-option -Wno-deprecated-register -Wno-header-guard -Wunused-const-variable -Wno-unused-function -Wno-pessimizing-move -Wno-missing-braces")

    if (NOT (${compiler_version} VERSION_LESS "3.3"))
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Woverloaded-shift-op-parentheses")
    endif ()

    if (NOT (${compiler_version} VERSION_LESS "3.4"))
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wheader-guard -Wlogical-not-parentheses -Wloop-analysis")
        #set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wunique-enum")
    endif ()

    # documentation of other possible warning flags from clang
    #set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Weverything -Wno-sign-conversion -Wno-weak-vtables -Wno-conversion -Wno-cast-align -Wno-padded -Wno-switch-enum -Wno-missing-noreturn -Wno-covered-switch-default -Wno-unreachable-code -Wno-global-constructors -Wno-exit-time-destructors")
endif()


# The NDEBUG macro is intentionally removed from release. One discussion on this is:
# http://www.drdobbs.com/an-exception-or-a-bug/184401686

# We want c++11, zlib and bzip2
set (CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS} ${CXX_WARN_FLAGS} -std=c++11")

if (GNU_COMPAT_COMPILER)
    set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set (CMAKE_CXX_FLAGS_RELEASE "-O3 -fomit-frame-pointer")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    #set (CMAKE_CXX_FLAGS_PROFILE "-O0 -g -pg -fprofile-arcs -ftest-coverage")
endif()


# add address sanitizer to debug mode:
set (USE_ADDRESS_SANITIZER false) # if true, turn on Address Sanitizer in debug for compilers which support this:

if (${USE_ADDRESS_SANITIZER})
    set (IS_ASAN false)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if (NOT (${compiler_version} VERSION_LESS "4.8"))
            set (IS_ASAN true)
        endif ()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        if (NOT (${compiler_version} VERSION_LESS "3.1"))
            set (IS_ASAN true)
        endif ()
    endif ()

    if (${IS_ASAN})
        set (CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-omit-frame-pointer ${CMAKE_CXX_FLAGS_DEBUG}")
    endif ()
endif ()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")

if (GNU_COMPAT_COMPILER)
  if (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[67]86$")
    ##
    ## Use scalar floating point instructions from the SSE instruction set.
    ## Note: Pentium3 SSE supports only single precision arithmetics
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -mfpmath=sse")
  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[345]86$")
    ##
    ## Prevent using 80bits registers (more consistent rounding)
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
  endif ()

endif()

