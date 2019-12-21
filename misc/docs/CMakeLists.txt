cmake_minimum_required(VERSION 3.13)

project(libqif)

find_package(Sphinx)

if(NOT SPHINX_EXECUTABLE)
	message(STATUS "Sphinx not installed, skipping docs ${SPHINX_EXECUTABLE}")
	add_custom_target(docs COMMAND echo Sphinx not installed)
	return()
endif()

include(../cmake/UseSphinxDoc.cmake)

add_sphinx_doc(
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	PYTHONPATH ${CMAKE_CURRENT_BINARY_DIR}/../../python_lib
	TARGET_NAME docs
)
add_dependencies(docs qif_module)
