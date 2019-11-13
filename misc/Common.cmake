include_guard(GLOBAL)

# Common configuration, included both in the library build, and in a standalone samples build.

# default type is Release, change with cmake -DCMAKE_BUILD_TYPE=
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# use unity builds by default for all targets EXCEPT osqp above (which breaks under unity builds)
#
if(NOT CMAKE_UNITY_BUILD AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
	set(CMAKE_UNITY_BUILD ON)
	set(DEFAULT_UNITY_BUILD ON)
	# set(CMAKE_UNITY_BUILD_BATCH_SIZE 0)
endif()

# use c++17
set(CMAKE_CXX_STANDARD 17)												# easier method for setting c++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
if(APPLE)
	LINK_DIRECTORIES(/usr/local/lib)
endif()

# which libraries to link
set(QIF_LIBS gmp gmpxx gsl gslcblas armadillo)

# use ortools, if available
#
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")	# get cmake modules from misc
find_package(ortools)
if(ortools_FOUND)
	message(STATUS "Found ortools")
	set(QIF_USE_ORTOOLS 1)									# this will be used in qif_bits/config.h
	list(APPEND QIF_LIBS ortools)							# add ortools to the list of linked libs
endif()

# use glpk, if available
find_library(LIB_GLPK NAME glpk PATH_SUFFIXES lib/)
if(NOT ${LIB_GLPK} STREQUAL "LIB_GLPK-NOTFOUND")
	message(STATUS "Found glpk")
	set(QIF_USE_GLPK 1)										# this will be used in qif_bits/config.h
	list(APPEND QIF_LIBS ${LIB_GLPK})						# add glpk to the list of linked libs
endif()