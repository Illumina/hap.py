add_subdirectory(klib)
add_subdirectory(jsoncpp)

set (KLIB_LIBRARY klib PARENT_SCOPE)
set (JSONCPP_LIBRARY jsoncpp PARENT_SCOPE)

install (DIRECTORY "htslib/htslib/" DESTINATION include/htslib FILES_MATCHING PATTERN "*.h")
install (DIRECTORY "intervaltree" DESTINATION include FILES_MATCHING PATTERN "*.h")
