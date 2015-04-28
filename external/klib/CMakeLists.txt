
# exclude things with external dependencies, like kurl
# file(GLOB LIBRARY_SOURCES *.c)
set(LIBRARY_SOURCES ksw.c)

set_source_files_properties(ksw.c PROPERTIES COMPILE_FLAGS -msse2)

add_library(klib STATIC ${LIBRARY_SOURCES})

set (KLIB_LIBRARY klib PARENT_SCOPE)

install (FILES "ksw.h" DESTINATION include)