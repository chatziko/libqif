cmake_minimum_required(VERSION 2.6)

project(libqif)


add_subdirectory(../external/gtest-1.7.0 gtest)						# build gtest in dir 'gtest'

include_directories(../external/gtest-1.7.0/include)

file(GLOB TEST_SOURCES *.cpp)										# get all *.cpp files in TEST_SOURCES

add_executable(run ${TEST_SOURCES})
target_link_libraries(run qif gtest_main)							# crete 'run' executable, link against libqif and gtest (with main function)

add_precompiled_header(run ../include/precompiled_tests.h)


# for 'make tests'
add_custom_target(tests)
add_dependencies(tests run)
