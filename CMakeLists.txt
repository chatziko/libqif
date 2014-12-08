cmake_minimum_required(VERSION 2.6)

project(libqif)

include(misc/PrecompiledHeader.cmake)									# for add_precompiled_header

add_definitions(-std=c++11 -Wall) #-fno-elide-constructor -O3 -march=natives

include_directories(include)


# build shard lib
#
file(GLOB LIB_SOURCES src/*.cpp)										# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})									# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk gmp gmpxx armadillo)						# libraries to link

add_precompiled_header(qif include/precompiled.h)


# other targets
#
add_subdirectory(tests EXCLUDE_FROM_ALL)								# build tests

add_subdirectory(samples EXCLUDE_FROM_ALL)								# build samples

add_subdirectory(misc/doxygen EXCLUDE_FROM_ALL)							# build doc

add_custom_target(allcode)												# 'all' builds just the library, 'allcode' builds also tests and samples
add_dependencies(allcode all tests samples)



# Print all variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

