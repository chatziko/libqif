cmake_minimum_required(VERSION 2.6)

project(libqif)

add_subdirectory(../external/backward-cpp backward-cpp)

MACRO(SUBDIRLIST result curdir)
	FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
	SET(dirlist "")
	FOREACH(child ${children})
		IF(IS_DIRECTORY ${curdir}/${child})
			LIST(APPEND dirlist ${child})
		ENDIF()
	ENDFOREACH()
	SET(${result} ${dirlist})
ENDMACRO()


add_custom_target(samples)												# "samples" target, depends on all samples

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach(DIR ${SUBDIRS})
	# for every dir 'foo' in samples, create a 'foo' executable by compiling samples/foo/*.cpp
	#
	file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/*.cpp)		# get all *.cpp files in SUBDIR
	set(DEST ${DIR})

	# make aux/<dest>-merged.cpp containing #include's for all source files
	# much faster than independent .cpp compilation, it avoids multiple template instantiation
	#
	file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/aux/)
	set(MERGED ${PROJECT_BINARY_DIR}/aux/${DEST}-merged.cpp)
	add_custom_command(
		OUTPUT ${MERGED}
		COMMAND perl -e 'print qq{\#include \"$$_\" \\n} for @ARGV' ${SOURCES} > ${MERGED}
		DEPENDS ${SOURCES}
	)

	add_executable(${DEST} ${MERGED} ${BACKWARD_ENABLE})
	target_link_libraries(${DEST} qif)									# link sample against libqif
	add_backward(${DEST})												# add stuff required by backward-cpp

	add_precompiled_header(${DEST} ../include/precompiled.h)

	add_dependencies(samples ${DEST})
endforeach()
