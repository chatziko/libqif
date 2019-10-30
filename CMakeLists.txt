cmake_minimum_required(VERSION 3.1)
cmake_policy(VERSION 3.1)												# avoid MACOSX_RPATH warning

project(libqif)

# default type is Release, change with cmake -DCMAKE_BUILD_TYPE=
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# build osqp in dir 'osqp', don't include targets in 'all' (we'll only depend on osqpstatic)
add_subdirectory(external/osqp osqp EXCLUDE_FROM_ALL)

# use unity builds by default for all targets EXCEPT osqp above (which breaks under unity builds)
#
if(NOT CMAKE_UNITY_BUILD AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
	set(CMAKE_UNITY_BUILD ON)
	# set(CMAKE_UNITY_BUILD_BATCH_SIZE 0)
endif()

# use c++17
set(CMAKE_CXX_STANDARD 17)												# easier method for setting c++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
if(APPLE)
	LINK_DIRECTORIES(/usr/local/lib)
endif()

# build shared lib
#
file(GLOB LIB_SOURCES src/*.cpp)											# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})										# libqif, depends on all src/*.cpp files
target_link_libraries(qif glpk gmp gmpxx gsl gslcblas armadillo osqpstatic)	# libraries to link

# enable precompiled headers in cmake >= 3.16.0
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
	target_precompile_headers(qif PRIVATE include/precompiled.h)
endif()

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
# NOTE: -fPIC and qif_EXPORTS are set automatically by cmake. We explicitely list them here as PUBLIC, so that they
#       are inherited by all targets that depend on qif, so that qif's precompiled header can be reused for them.
#
target_compile_options(qif PUBLIC -Wall -Wextra -pedantic  $<$<CONFIG:Release>:-march=native> -fPIC)
target_compile_definitions(qif PUBLIC $<$<CONFIG:Release>:ARMA_NO_DEBUG> qif_EXPORTS)

# includes and config file
configure_file(include/qif_bits/config.h include/qif_bits/config.h)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/include" include)

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

install(TARGETS qif
	ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
	LIBRARY DESTINATION ${INSTALL_LIB_DIR})

# Print all variables
if(VERBOSE)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endif()

