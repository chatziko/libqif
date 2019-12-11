cmake_minimum_required(VERSION 3.13)

project(libqif)

include(misc/Common.cmake)

# build osqp in dir 'osqp', don't include targets in 'all' (we'll only depend on osqpstatic)
set(CMAKE_UNITY_BUILD OFF)								# never use unity for osqp, it breaks
add_subdirectory(external/osqp osqp EXCLUDE_FROM_ALL)
set(CMAKE_UNITY_BUILD ${DEFAULT_UNITY_BUILD})

# build mp++ in dir 'mppp'
add_subdirectory(external/mppp mppp EXCLUDE_FROM_ALL)

# build shared lib
#
file(GLOB LIB_SOURCES src/*.cpp)						# get all src/*.cpp files in LIB_SOURCES
add_library(qif SHARED ${LIB_SOURCES})					# libqif, depends on all src/*.cpp files
target_link_libraries(qif osqpstatic ${QIF_LIBS})		# libraries to link (QIF_LIBS is set in Common.cmake)

# enable precompiled headers in cmake >= 3.16.0
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
	target_precompile_headers(qif PRIVATE include/precompiled.h)
endif()

# compiler options
# NOTE: -fPIC and qif_EXPORTS are set automatically by cmake. We explicitely list them here as PUBLIC, so that they
#       are inherited by all targets that depend on qif, so that qif's precompiled header can be reused for them.
#
target_compile_options(qif PUBLIC -Wall -Wextra -pedantic  $<$<CONFIG:Release>:-march=native> -fPIC -pthread)
target_compile_definitions(qif PUBLIC $<$<CONFIG:Release>:ARMA_NO_DEBUG> qif_EXPORTS)

# includes and config file
configure_file(include/qif_bits/config.h include/qif_bits/config.h)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/include" include)

# other targets
#
add_subdirectory(samples EXCLUDE_FROM_ALL)								# build samples

add_subdirectory(python EXCLUDE_FROM_ALL)								# build python bindings

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
	PATTERN "miniball" EXCLUDE
	PATTERN "_*" EXCLUDE)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION "include")

# install mp++'s includes, they have to be globally installed cause mp++ uses #include<...>
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/mppp/include/ DESTINATION "include")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mppp/include/ DESTINATION "include")

# install OSQP's includes under qif_bits/osqp
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/osqp/include/ DESTINATION "include/qif_bits/osqp")

# install miniball's includes under qif_bits/miniball
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/miniball/cpp/main/ DESTINATION "include/qif_bits/miniball")

install(TARGETS qif mp++
	ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
	LIBRARY DESTINATION ${INSTALL_LIB_DIR})

# Print all variables
if(VERBOSE)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endif()

