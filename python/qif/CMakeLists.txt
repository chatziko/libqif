cmake_minimum_required(VERSION 2.6)

project(libqif)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(FILE ${SOURCES})
	get_filename_component(BASE ${FILE} NAME_WE) 

	pybind11_add_module(${BASE} ${FILE})

	target_link_libraries(${BASE} PRIVATE qif)

	# enable precompiled headers in cmake >= 3.16.0, EXCEPT on linux/clang
	if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0" AND NOT (UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
		target_precompile_headers(${BASE} REUSE_FROM qif)
	endif()

	add_dependencies(python ${BASE})
endforeach()
