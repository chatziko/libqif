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


# dummy target, used only to generate a common precompiled header for all samples
file(WRITE "${PROJECT_BINARY_DIR}/empty.cpp" "int main(){}")
add_executable(samples_pch "${PROJECT_BINARY_DIR}/empty.cpp")
set_target_properties(samples_pch PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "../include/precompiled.h")
cotire(samples_pch)


add_custom_target(samples)							# "samples" target, depends on all samples

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach(DIR ${SUBDIRS})
	# for every dir 'foo' in samples, create a 'foo' executable by compiling samples/foo/*.cpp
	#
	file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/*.cpp)		# get all *.cpp files in SUBDIR
	set(DEST ${DIR})

	add_executable(${DEST} ${SOURCES} ${BACKWARD_ENABLE})
	target_link_libraries(${DEST} qif)									# link sample against libqif
	add_backward(${DEST})												# add stuff required by backward-cpp

	# we manually add the dependency and the -include flag to the common precompiled header
	cotire(${DEST})
	add_dependencies("${DEST}_unity" samples_pch)
	target_compile_options("${DEST}_unity" PRIVATE -include ${PROJECT_BINARY_DIR}/cotire/samples_pch_CXX_prefix.hxx)

	add_dependencies(samples "${DEST}_unity")
endforeach()
