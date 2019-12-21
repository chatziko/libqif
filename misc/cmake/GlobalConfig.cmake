include_guard(GLOBAL)

# Global config, included in all (potentially) top-level CMakeLists: qif / samples / python

# default type is Release, change with cmake -DCMAKE_BUILD_TYPE=
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# use unity builds by default for all targets
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
	link_directories(/usr/local/lib)
	include_directories(/usr/local/include)
endif()

# without this the installed library has empty rpath, which confuses delocate
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# get QIF_VERSION from the most recent vN.N.N tag
execute_process(COMMAND git describe --tags OUTPUT_VARIABLE GIT_TAG)
if(${GIT_TAG} MATCHES "v([0-9.]+)")
	set(QIF_VERSION ${CMAKE_MATCH_1})
else()
	set(QIF_VERSION "dev")
endif()
message(STATUS "Building QIF Version: ${QIF_VERSION}")

# find ortools
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")	# get cmake modules from misc
find_package(ortools)
if(ortools_FOUND)
	message(STATUS "Found ortools")
	set(QIF_USE_ORTOOLS 1)									# this will be used in qif_bits/config.h
endif()

# find glpk
find_library(LIB_GLPK NAME glpk PATH_SUFFIXES lib/)
if(NOT ${LIB_GLPK} STREQUAL "LIB_GLPK-NOTFOUND")
	message(STATUS "Found glpk: ${LIB_GLPK}")
	set(QIF_USE_GLPK 1)										# this will be used in qif_bits/config.h
endif()

# Macros
MACRO(SUBDIRLIST result curdir)
	FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
	SET(dirlist "")
	FOREACH(child ${children})
		IF(IS_DIRECTORY ${curdir}/${child})
			LIST(APPEND dirlist ${child})
		ENDIF()
	ENDFOREACH()
	SET(${result} ${dirlist})
ENDMACRO()