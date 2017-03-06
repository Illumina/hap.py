#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
# 
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
# - Try to find htslib
# Define the environment variable HTSLIB_ROOT to point to a htslib 
# download directory.
# Once done this will define
#  HTSLIB_FOUND - System has htslib
#  HTSLIB_INCLUDE_DIRS - The htslib include directories
#  HTSLIB_LIBRARIES - The libraries needed

# message (STATUS "FindHTSlib.cmake")
# message (STATUS "ENV HTSLIB_ROOT : " $ENV{HTSLIB_ROOT})

if( DEFINED ENV{HTSLIB_ROOT} )
    set( HTSLIB_ROOT $ENV{HTSLIB_ROOT} )
    # message (STATUS "HTSLIB_ROOT : ${HTSLIB_ROOT}")

    # set library directories
    set(HTSLIB_LIBRARY_DIR ${HTSLIB_ROOT})
    # set include directories
    set(HTSLIB_INCLUDE_DIR ${HTSLIB_ROOT})

    find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
              HINTS ${HTSLIB_INCLUDEDIR}
              )

  	find_library(HTSLIB_LIBRARY NAMES libhts.a
                 HINTS ${HTSLIB_LIBRARY_DIR})

    set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY} )
    set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR} )

    set(HTSLIB_FOUND TRUE)
else()
    find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
              HINTS "${CMAKE_BINARY_DIR}/include")

	  find_library(HTSLIB_LIBRARY NAMES libhts.a
                 HINTS "${CMAKE_BINARY_DIR}/lib")

    set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY} )
    set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR} )

    set(HTSLIB_FOUND TRUE)
endif()

include_directories(${HTSLIB_INCLUDE_DIRS})

message (STATUS "HTSLIB_INCLUDE_DIR=" ${HTSLIB_INCLUDE_DIR})
#message (STATUS "HTSLIB_LIBRARIES=" ${HTSLIB_LIBRARIES})
