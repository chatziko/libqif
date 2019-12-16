include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.13)

project(libqif)

include(misc/cmake/GlobalConfig.cmake)						# global cmake config

# various targets
#
if(NOT DEFINED BUILD_QIF)
	set(BUILD_QIF true)										# build by default
endif()

add_subdirectory(lib)										# build library
add_subdirectory(samples		EXCLUDE_FROM_ALL)			# build samples
add_subdirectory(python_lib		EXCLUDE_FROM_ALL)			# build python lib
add_subdirectory(tests			EXCLUDE_FROM_ALL)			# build tests, NEEDS TO BE AFTER python_lib, cause it somehow messes with python detection!
add_subdirectory(misc/doxygen	EXCLUDE_FROM_ALL)			# build doc

add_custom_target(allcode DEPENDS samples qif_module tests)	# 'all' builds just the library, 'allcode' builds also tests, samples and python bindings


# for debugging, print all variables
if(VERBOSE)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endif()
