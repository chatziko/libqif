cmake_minimum_required(VERSION 2.6)

project(libqif)

include(../python/cmake/GlobalConfig.cmake)
include(../python/cmake/QifConfig.cmake)

# use backward.cpp of pretty stacktraces
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../external/backward-cpp/CMakeLists.txt")
	set(CMAKE_UNITY_BUILD OFF)									# never use unity build for backward.cpp, it breaks
	add_subdirectory(../external/backward-cpp backward-cpp)
	set(CMAKE_UNITY_BUILD ${DEFAULT_UNITY_BUILD})

else()
	# backward.cpp not present, add dummy commands to avoid errors
	set(BACKWARD_ENABLE "")
	macro(add_backward target)
	endmacro()
endif()

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
	if(NOT SYSTEM_QIF AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
		target_precompile_headers(${DEST} REUSE_FROM qif)
	endif()

	add_dependencies(samples ${DEST})
endforeach()
