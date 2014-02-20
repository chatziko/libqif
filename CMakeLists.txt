cmake_minimum_required(VERSION 2.6)

project(libqif)

add_definitions(-std=c++11)

include_directories(include)

# build shard lib
#
file(GLOB LIB_SOURCES src/*.cpp)										# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})									# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk)											# link libqif against libglpk


add_subdirectory(tests EXCLUDE_FROM_ALL)								# build tests

add_subdirectory(samples EXCLUDE_FROM_ALL)								# build samples

add_subdirectory(doxygen EXCLUDE_FROM_ALL)								# build doc

add_custom_target(allcode)												# 'all' builds just the library, 'allcode' builds also tests and samples
add_dependencies(allcode all tests samples)

