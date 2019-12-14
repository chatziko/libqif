include_guard(GLOBAL)

# qif library config, included in all (potentially) top-level CMakeLists: qif / samples / python


# create qif target
if(${BUILD_QIF})
	file(GLOB LIB_SOURCES src/*.cpp)							# get all src/*.cpp files in LIB_SOURCES
	add_library(qif SHARED ${LIB_SOURCES})						# libqif, depends on all src/*.cpp files

	target_include_directories(qif PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include" include)

	# enable precompiled headers in cmake >= 3.16.0
	if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
		target_precompile_headers(qif PRIVATE include/precompiled.h)
	endif()

else()
	# using existing qif
	find_library(SYSTEM_QIF NAME qif PATH_SUFFIXES lib/)
	add_library(qif INTERFACE)
	target_link_libraries(qif INTERFACE ${SYSTEM_QIF})			# targets using qif will link to the system qif
endif()

# linked libraries
#
# PRIVATE, only needed for qif itself
if(${BUILD_QIF})
	target_link_libraries(qif PRIVATE osqpstatic gsl gslcblas mp++)		# mp++ is not really needed, we only list it to get its include_directories
endif()

# INTERFACE, needed for targets _using_ qif
target_link_libraries(qif INTERFACE gmp mp++ armadillo)
# if(UNIX AND NOT APPLE)
	# dl and rt are needed in manylinux
	#target_link_libraries(qif INTERFACE ${CMAKE_DL_LIBS} rt)
# endif()

# compiler options
#
# NOTE: -fPIC and qif_EXPORTS are set automatically by cmake. We explicitely list them here as INTERFACE, so that they
#       are inherited by all targets that depend on qif, so that qif's precompiled header can be reused for them.
#
target_compile_options(qif INTERFACE -Wall -Wextra -pedantic  $<$<CONFIG:Release>:-march=native> -fPIC -pthread)
target_compile_definitions(qif INTERFACE $<$<CONFIG:Release>:ARMA_NO_DEBUG> qif_EXPORTS)

if(${BUILD_QIF})
	# if we're building qif, copy to PRIVATE so that they are used also for qif itself
	get_target_property(TMP1 qif INTERFACE_COMPILE_DEFINITIONS)
	get_target_property(TMP2 qif INTERFACE_COMPILE_OPTIONS)
	target_compile_definitions(qif PRIVATE ${TMP1})
	target_compile_options    (qif PRIVATE ${TMP2})
endif()


# use ortools, if available
#
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")	# get cmake modules from misc
find_package(ortools)
if(ortools_FOUND)
	message(STATUS "Found ortools")
	set(QIF_USE_ORTOOLS 1)									# this will be used in qif_bits/config.h
	target_link_libraries(qif INTERFACE ortools)			# add ortools to the list of linked libs
endif()

# use glpk, if available
find_library(LIB_GLPK NAME glpk PATH_SUFFIXES lib/)
if(NOT ${LIB_GLPK} STREQUAL "LIB_GLPK-NOTFOUND")
	message(STATUS "Found glpk")
	set(QIF_USE_GLPK 1)										# this will be used in qif_bits/config.h
	target_link_libraries(qif INTERFACE ${LIB_GLPK})		# add glpk to the list of linked libs
endif()
