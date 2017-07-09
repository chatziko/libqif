cmake_minimum_required(VERSION 2.6)

project(libqif)


add_subdirectory(../external/googletest googletest)					# build googletest in dir 'googletest'

include_directories(../external/googletest/include)

file(GLOB TEST_SOURCES *.cpp)										# get all *.cpp files in TEST_SOURCES

add_executable(run ${TEST_SOURCES})
target_link_libraries(run qif gtest)								# crete 'run' executable, link against libqif and googletest

set_target_properties(run PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "../include/precompiled_tests.h")
cotire(run)
add_dependencies(run_unity run_pch)									# cotire seems to add the -include <precompiled> flag in the _unity build, but not a dependency to the _pch target that builds <precompiled>


# for 'make tests'
add_custom_target(tests)
add_dependencies(tests run_unity)
