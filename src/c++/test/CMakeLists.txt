file(GLOB TEST_SOURCES test_*.cpp test_*.c)

add_executable(test_haplotypes ${TEST_SOURCES})
target_link_libraries(test_haplotypes ${HAPLOTYPES_ALL_LIBS})
SET(LINK_LIBRARIES test_haplotypes ${Boost_TEST_EXEC_MONITOR_LIBRARY})
