set(LIBRARY_SOURCES jsoncpp.cpp)

# no warnings. not ideal, but -Wno-enum-compare doesn't seem to do the trick
set_source_files_properties(jsoncpp.cpp PROPERTIES COMPILE_FLAGS "-w")

add_library(jsoncpp STATIC ${LIBRARY_SOURCES})

set (JSONCPP_LIBRARY jsoncpp PARENT_SCOPE)

install (DIRECTORY "json/" DESTINATION include/json FILES_MATCHING PATTERN "*.h")