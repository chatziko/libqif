cmake_minimum_required(VERSION 2.6)

project(libqif)

include(misc/PrecompiledHeader.cmake)									# for add_precompiled_header

# default type is Release, change with cmake -DCMAKE_BUILD_TYPE=
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# compiler flags
#
add_definitions(-std=c++11 -Wall -Wextra -pedantic)						# common flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	message(STATUS "Clang compiler on MacOS X detected. Added '-stdlib=libc++' to compiler flags")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -DARMA_NO_DEBUG")		# Release flags

include_directories(include)


# build shard lib
#
file(GLOB LIB_SOURCES src/*.cpp)										# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})									# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk gmp gmpxx gsl armadillo)					# libraries to link

add_precompiled_header(qif include/precompiled.h)


# other targets
#
add_subdirectory(samples EXCLUDE_FROM_ALL)								# build samples

add_subdirectory(misc/doxygen EXCLUDE_FROM_ALL)							# build doc

add_custom_target(allcode)												# 'all' builds just the library, 'allcode' builds also tests and samples
add_dependencies(allcode samples)

# tests only external/googletest is checked out
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/CMakeLists.txt")
	add_subdirectory(tests EXCLUDE_FROM_ALL)							# build tests
	add_dependencies(allcode tests)
endif()



# Print all variables
if(VERBOSE)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endif()

