include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.13)

project(libqif)

# create qif target that either builds shared lib, or uses the system one
#
if("${BUILD_QIF}")
	# build osqp in dir 'osqp', don't include targets in 'all' (we'll only depend on osqpstatic)
	set(CMAKE_UNITY_BUILD OFF)								# never use unity for osqp, it breaks
	add_subdirectory(../external/osqp osqp EXCLUDE_FROM_ALL)
	set(CMAKE_UNITY_BUILD ${DEFAULT_UNITY_BUILD})

	# build mp++ in dir 'mppp'
	add_subdirectory(../external/mppp mppp EXCLUDE_FROM_ALL)

	file(GLOB LIB_SOURCES src/*.cpp)							# get all src/*.cpp files in LIB_SOURCES
	add_library(qif SHARED ${LIB_SOURCES})						# libqif, depends on all src/*.cpp files

	target_include_directories(qif PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include" include)

	# enable precompiled headers in cmake >= 3.16.0
	if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
		target_precompile_headers(qif PRIVATE include/precompiled.h)
	endif()

else()
	# using existing qif
	find_library(QIF_PATH NAME qif PATH_SUFFIXES lib/)
	if(${QIF_PATH} STREQUAL "QIF_PATH-NOTFOUND")
		message(FATAL_ERROR "QIF_BUILD is not set, but libqif is not installed in the system")
	endif()

	add_library(qif INTERFACE)
	target_link_libraries(qif INTERFACE ${QIF_PATH})					# targets using qif will link to the system qif
endif()

# linked libraries
#
# PRIVATE, only needed for qif itself
if("${BUILD_QIF}")
	target_link_libraries(qif PRIVATE osqpstatic gsl gslcblas mp++)		# mp++ is not really needed, we only list it to get its include_directories
endif()

# INTERFACE, needed for targets _using_ qif
target_link_libraries(qif INTERFACE gmp mp++ armadillo)
if(${QIF_USE_ORTOOLS})
	target_link_libraries(qif INTERFACE ortools)
endif()
if(${QIF_USE_GLPK})
	target_link_libraries(qif INTERFACE ${LIB_GLPK})
endif()

# compiler options
#
# NOTE: -fPIC and qif_EXPORTS are set automatically by cmake. We explicitely list them here as INTERFACE, so that they
#       are inherited by all targets that depend on qif, so that qif's precompiled header can be reused for them.
#
target_compile_options(qif INTERFACE -Wall -Wextra -pedantic  $<$<CONFIG:Release>:-march=native> -fPIC -pthread)
target_compile_definitions(qif INTERFACE $<$<CONFIG:Release>:ARMA_NO_DEBUG> qif_EXPORTS)

if("${BUILD_QIF}")
	# if we're building qif, copy to PRIVATE so that they are used also for qif itself
	get_target_property(TMP1 qif INTERFACE_COMPILE_DEFINITIONS)
	get_target_property(TMP2 qif INTERFACE_COMPILE_OPTIONS)
	target_compile_definitions(qif PRIVATE ${TMP1})
	target_compile_options    (qif PRIVATE ${TMP2})
endif()

# config file
configure_file(include/qif_bits/config.h include/qif_bits/config.h)



# installation
#
if("${BUILD_QIF}")	# can only install if we built
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
		PATTERN "precompiled.h" EXCLUDE
		PATTERN "osqp" EXCLUDE
		PATTERN "miniball" EXCLUDE
		PATTERN "_*" EXCLUDE)
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION "include")

	# install mp++'s includes, they have to be globally installed cause mp++ uses #include<...>
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../external/mppp/include/ DESTINATION "include")
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mppp/include/ DESTINATION "include")

	# install OSQP's includes under qif_bits/osqp
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../external/osqp/include/ DESTINATION "include/qif_bits/osqp")

	# install miniball's includes under qif_bits/miniball
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../external/miniball/cpp/main/ DESTINATION "include/qif_bits/miniball")

	install(TARGETS qif mp++
		ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
		LIBRARY DESTINATION ${INSTALL_LIB_DIR})
endif()

