cmake_minimum_required(VERSION 3.1)
cmake_policy(VERSION 3.1)												# avoid MACOSX_RPATH warning
if(POLICY CMP0058)
	cmake_policy(SET CMP0058 NEW)										# void cotire warning under ninja
endif()

project(libqif)

include(misc/cotire.cmake)												# for unity builds and precompiled headers
set (COTIRE_MINIMUM_NUMBER_OF_TARGET_SOURCES 1)

# default type is Release, change with cmake -DCMAKE_BUILD_TYPE=
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# build osqp in dir 'osqp', don't include targets in 'all' (we'll only depend on osqpstatic)
add_subdirectory(external/osqp osqp EXCLUDE_FROM_ALL)

# use c++11
set(CMAKE_CXX_STANDARD 11)											# easier method for setting c++11, also avoid cotire issue with ortools
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# build shared lib
#
file(GLOB LIB_SOURCES src/*.cpp)										# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED EXCLUDE_FROM_ALL ${LIB_SOURCES})					# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk gmp gmpxx gsl gslcblas armadillo osqpstatic)		# libraries to link

# use ortools, if available
#
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/misc")			# get cmake modules from misc
find_package(ortools)
if(ortools_FOUND)
	message(STATUS "Found ortools")
	set(QIF_USE_ORTOOLS 1)							# this will be used in qif_bits/config.h
	target_link_libraries(qif ortools::ortools)
endif()

# compiler options
#
target_compile_options(qif PUBLIC -Wall -Wextra -pedantic  $<$<CONFIG:Release>:-march=native>)
target_compile_definitions(qif PUBLIC $<$<CONFIG:Release>:ARMA_NO_DEBUG>)

# includes and config file
configure_file(include/qif_bits/config.h include/qif_bits/config.h)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/include" include)

# apply cotire
#
set_target_properties(qif PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "include/precompiled.h")
cotire(qif)
add_dependencies(qif_unity qif_pch)										# cotire seems to add the -include <precompiled> flag in the _unity build, but not a dependency to the _pch target that builds <precompiled>
set_target_properties(qif_unity PROPERTIES EXCLUDE_FROM_ALL FALSE)		# 'all' builds qif_unity

# other targets
#
add_subdirectory(samples EXCLUDE_FROM_ALL)								# build samples

add_subdirectory(misc/doxygen EXCLUDE_FROM_ALL)							# build doc

add_custom_target(allcode DEPENDS samples)								# 'all' builds just the library, 'allcode' builds also tests and samples

# tests, only if external/googletest is checked out
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/CMakeLists.txt")
	add_subdirectory(tests EXCLUDE_FROM_ALL)							# build tests
	add_dependencies(allcode tests)
endif()

# installation
#
if(NOT INSTALL_LIB_DIR) # Allow for the "lib" directory to be specified on the command line
	set(INSTALL_LIB_DIR "lib")
	if(UNIX AND NOT APPLE)	 # I don't know how Mac OS handles 64 bit systems
		if(CMAKE_SIZEOF_VOID_P EQUAL 8)
			message(STATUS "*** Detected 64 bit system")
			# use "lib64" only on systems that have it (eg. Fedora, Red Hat).
			# if installing in /usr/local/, the use of lib64 might be unreliable on systems which have /usr/local/lib64 but don't use it
			if(IS_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib64")
				unset(INSTALL_LIB_DIR)
				set(INSTALL_LIB_DIR "lib64")
				message(STATUS "*** ${CMAKE_INSTALL_PREFIX}/lib64/ exists, so destination directory for the run-time library changed to ${CMAKE_INSTALL_PREFIX}/lib64/")
				message(STATUS "*** Your system and/or compiler must search ${CMAKE_INSTALL_PREFIX}/lib64/ during linking")
			endif()
		endif()
	endif()
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION "include"
	PATTERN "precompiled*.h" EXCLUDE
	PATTERN "tests_aux.h" EXCLUDE
	PATTERN "osqp" EXCLUDE
	PATTERN "_*" EXCLUDE)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION "include")

# install OSQP's includes under qif_bits/osqp
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/osqp/include/ DESTINATION "include/qif_bits/osqp")

install(TARGETS qif_unity
	ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
	LIBRARY DESTINATION ${INSTALL_LIB_DIR})

# Print all variables
if(VERBOSE)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endif()

