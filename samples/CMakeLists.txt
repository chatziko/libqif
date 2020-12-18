include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.13)

include(../misc/cmake/GlobalHead.cmake)
project(libqif)
include(../misc/cmake/GlobalConfig.cmake)
include(../misc/cmake/BackwardCpp.cmake)

# If this dir is the root, we need to include qif (if not qif is already loaded, and add_subdirectory will do nothing)
add_subdirectory(../lib lib)

add_custom_target(samples)							# "samples" target, depends on all samples

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach(DIR ${SUBDIRS})
	# for every dir 'foo' in samples, create a 'foo' executable by compiling samples/foo/*.cpp
	#
	file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/*.cpp)		# get all *.cpp files in SUBDIR
	set(DEST ${DIR})

	add_executable(${DEST} ${SOURCES} ${BACKWARD_ENABLE})
	target_link_libraries(${DEST} qif)								# link sample against libqif
	add_backward(${DEST})											# add stuff required by backward-cpp

	# enable precompiled headers in cmake >= 3.16.0
	if("${BUILD_QIF}" AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
		target_precompile_headers(${DEST} REUSE_FROM qif)
	endif()

	add_dependencies(samples ${DEST})
endforeach()
