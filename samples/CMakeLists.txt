cmake_minimum_required(VERSION 2.6)

project(libqif)

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

	add_executable(${DEST} ${SOURCES})
	target_link_libraries(${DEST} qif)									# link sample against libqif

	add_precompiled_header(${DEST} ../include/precompiled.h)

	add_dependencies(samples ${DEST})
endforeach()

