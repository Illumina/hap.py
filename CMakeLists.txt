# Top-level Cmake file for Haplotype Tools
#
# Author: Peter Krusche <pkrusche@illumina.com>
#
cmake_minimum_required (VERSION 2.8)
project (HAPLOTYPES)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(cxx)
include(configureFiles)

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)

# un-break library finding
set(Boost_NO_BOOST_CMAKE         ON)
set(Boost_NO_SYSTEM_PATHS        ON)

find_package(Boost 1.50.0 COMPONENTS thread iostreams regex unit_test_framework filesystem system program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories (${HAPLOTYPES_SOURCE_DIR}/external/zlib-1.2.8)
set(ZLIB_LIBRARIES ${HAPLOTYPES_SOURCE_DIR}/external/zlib-1.2.8/libz.a)

include_directories (${HAPLOTYPES_SOURCE_DIR}/external/klib)
include_directories (${HAPLOTYPES_SOURCE_DIR}/external/intervaltree)
include_directories (${HAPLOTYPES_SOURCE_DIR}/external/jsoncpp)
include_directories (${HAPLOTYPES_SOURCE_DIR}/src/c++/include)

# install bcftools, samtools, htslib

option(USE_MUSCLE "Build Muscle3 Multiple Aligner (required for approximate haplotype matching)" OFF)
option(USE_SGE    "Require the --force-interactive option to run interactively in hap.py ()" OFF)

if(USE_MUSCLE)
    message("Muscle3/approximate haplotype matching is enabled")
    execute_process(
        COMMAND ./make_dependencies.sh build_muscle
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external)
    configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/muscle3.8.31/src" "muscle" "${CMAKE_BINARY_DIR}/bin")

    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAS_MUSCLE")
    set (MUSCLE_AVAILABLE 1)
else()
    message("Muscle3/approximate haplotype matching is disabled")
    execute_process(
        COMMAND ./make_dependencies.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external)
    configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/muscle3.8.31/src" "muscle" "${CMAKE_BINARY_DIR}/bin")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_MUSCLE")
    set (MUSCLE_AVAILABLE 0)
endif()

if(USE_SGE)
    message("SGE mode selected, you will need to use --force-interactive to run hap.py")
    set (SGE_AVAILABLE 1)
else()
    set (SGE_AVAILABLE 0)
endif()

configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/zlib-1.2.8" "libz*" "${CMAKE_BINARY_DIR}/lib")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/htslib" "libhts*" "${CMAKE_BINARY_DIR}/lib")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/htslib" "tabix" "${CMAKE_BINARY_DIR}/bin")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/htslib" "bgzip" "${CMAKE_BINARY_DIR}/bin")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/bcftools" "bcftools" "${CMAKE_BINARY_DIR}/bin")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}/external/samtools" "samtools" "${CMAKE_BINARY_DIR}/bin")

include(FindHTSLib)

find_package (Threads REQUIRED)

# make libraries first
add_subdirectory (external)
add_subdirectory (src/c++/lib)

set(HAPLOTYPES_ALL_LIBS ${HAPLOTYPES_LIBRARY}
                        ${GSSW_LIBRARY}
                        ${KLIB_LIBRARY}
                        ${JSONCPP_LIBRARY}
                        ${Boost_LIBRARIES}
                        ${HTSLIB_LIBRARIES}
                        ${ZLIB_LIBRARIES}
                        ${CMAKE_THREAD_LIBS_INIT})


include(GetGitRevisionDescription)
git_describe(HAPLOTYPES_VERSION)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/c++/include/Version.hh.in"
               "${CMAKE_BINARY_DIR}/include/Version.hh")

include_directories(${CMAKE_BINARY_DIR}/include)

add_subdirectory (src/c++/main)
add_subdirectory (src/data)
add_subdirectory (src/c++/test)
add_subdirectory (src/python)

# install


# dev install static libraries
install (DIRECTORY "${CMAKE_BINARY_DIR}/bin/" DESTINATION bin FILES_MATCHING PATTERN "*" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_READ WORLD_EXECUTE)
install (DIRECTORY "${CMAKE_BINARY_DIR}/lib/" DESTINATION lib)
install (DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/c++/include/" DESTINATION include)
