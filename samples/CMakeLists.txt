cmake_minimum_required(VERSION 2.6)

project(libqif)


add_custom_target(samples)												# "samples" target, depends on all samples

file(GLOB SAMPLE_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)	# get all *.cpp files, no path

foreach(SOURCE ${SAMPLE_SOURCES})
	string(REPLACE ".cpp" "" DEST ${SOURCE})							# SOURCE = foo.cpp, DEST = foo

	add_executable(${DEST} ${SOURCE})
	target_link_libraries(${DEST} qif)									# link sample against libqif

	add_dependencies(samples ${DEST})
endforeach()

