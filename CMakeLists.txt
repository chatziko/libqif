cmake_minimum_required(VERSION 2.6)

add_subdirectory(external/gtest-1.7.0)

project(libqif)

add_definitions(-std=c++11)

include_directories(inc)

# build shard lib
#
file(GLOB LIB_SOURCES src/*.cpp)				# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})			# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk)					# link libqif against libglpk


# build 'test' binary including all tests
#
include_directories(external/gtest-1.7.0/include)

file(GLOB TEST_SOURCES tests/*.cc)				# get all tests/*.cc files in TEST_SOURCES

add_executable(test ${TEST_SOURCES})
target_link_libraries(test qif gtest_main)		# crete 'test' executable, link against libqif and gtest (with main function)


# samples, one binary per .cpp files:
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "samples")

file(GLOB SAMPLE_SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/samples" "samples/*.cpp")	# *.cpp files in samples, without the path

foreach(SOURCE ${SAMPLE_SOURCES})
	string(REPLACE ".cpp" "" DEST ${SOURCE})	# SOURCE = foo.cpp, DEST = foo

	add_executable(${DEST} "samples/${SOURCE}")
	target_link_libraries(${DEST} qif)			# link sample against libqif
endforeach()

