cmake_minimum_required(VERSION 2.6)

add_subdirectory(external/gtest-1.7.0)

project(libqif)

add_definitions(-std=c++11)

include_directories(inc)

file(GLOB LIB_SOURCES src/*.cpp)				# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})			# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk)					# link libqif against libglpk

add_executable(sample samples/Main.cpp)			# sample prog
target_link_libraries(sample qif)				# link sample against libqif


# tests
#
include_directories(external/gtest-1.7.0/include)

file(GLOB TEST_SOURCES tests/*.cc)				# get all tests/*.cc files in TEST_SOURCES

add_executable(test ${TEST_SOURCES})
target_link_libraries(test qif gtest_main)		# crete 'test' executable, link against libqif and gtest (with main function)

