
file(GLOB_RECURSE LIBRARY_SOURCES *.cpp *.c)
add_library     (haplotypes STATIC ${LIBRARY_SOURCES})

set (HAPLOTYPES_LIBRARY haplotypes PARENT_SCOPE)