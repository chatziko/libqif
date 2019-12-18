cmake_minimum_required(VERSION 3.13)

project(libqif)

include(../misc/cmake/BackwardCpp.cmake)

add_custom_target(tests)												# for 'make tests'

if(EXISTS "../external/googletest/CMakeLists.txt")						# only if googletest is checked out

	add_subdirectory(../external/googletest googletest)					# build googletest in dir 'googletest'

	file(GLOB_RECURSE TEST_SOURCES *.cpp)								# get all *.cpp files in TEST_SOURCES

	add_executable(run ${TEST_SOURCES} ${BACKWARD_ENABLE})				# create 'run' executable
	target_link_libraries(run qif gtest)								# link against libqif and googletest
	target_include_directories(run PRIVATE include)
	add_backward(run)													# add stuff required by backward-cpp

	# enable precompiled headers in cmake >= 3.16.0
	if("${BUILD_QIF}" AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
		target_precompile_headers(run REUSE_FROM qif)
	endif()

	add_dependencies(tests run)

endif()
