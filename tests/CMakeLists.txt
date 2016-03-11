cmake_minimum_required(VERSION 2.6)

project(libqif)


add_subdirectory(../external/googletest googletest)					# build googletest in dir 'googletest'

include_directories(../external/googletest/include)

file(GLOB TEST_SOURCES *.cpp)										# get all *.cpp files in TEST_SOURCES

# make aux/run-merged.cpp containing #include's for all source files
# much faster than independent .cpp compilation, it avoids multiple template instantiation
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/aux/)
set(MERGED ${PROJECT_BINARY_DIR}/aux/run-merged.cpp)
add_custom_command(
	OUTPUT ${MERGED}
	COMMAND perl -e 'print qq{\#include \"$$_\" \\n} for @ARGV' ${TEST_SOURCES} > ${MERGED}
	DEPENDS ${TEST_SOURCES}
)

add_executable(run ${MERGED})
target_link_libraries(run qif gtest)								# crete 'run' executable, link against libqif and googletest

add_precompiled_header(run ../include/precompiled_tests.h)


# for 'make tests'
add_custom_target(tests)
add_dependencies(tests run)
