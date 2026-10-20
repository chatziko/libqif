include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.13)
enable_testing()

include(misc/cmake/GlobalHead.cmake)
project(libqif)
include(misc/cmake/GlobalConfig.cmake)						# global cmake config

# various targets
#
if(NOT DEFINED BUILD_QIF)
	set(BUILD_QIF true)										# build by default
endif()

add_subdirectory(lib_cpp)									# build cpp library

find_package(Python3 COMPONENTS Interpreter Development.Module)	# python is only needed for the python library
if(Python3_FOUND)
	add_subdirectory(lib_python	EXCLUDE_FROM_ALL)			# build python library
	set(PYTHON_TARGET qif_python)
else()
	message(STATUS "Python not found, skipping the python library")
endif()

add_subdirectory(samples		EXCLUDE_FROM_ALL)			# build samples
add_subdirectory(tests_cpp		EXCLUDE_FROM_ALL)			# build tests, NEEDS TO BE AFTER lib_python, cause it somehow messes with python detection!
add_subdirectory(misc/docs		EXCLUDE_FROM_ALL)			# build doc

add_custom_target(allcode DEPENDS samples ${PYTHON_TARGET} tests_cpp)	# 'all' builds just the library, 'allcode' builds also tests, samples and python bindings


# for debugging, print all variables
if(VERBOSE)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endif()
