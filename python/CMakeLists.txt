cmake_minimum_required(VERSION 2.6)

project(libqif)

include(../python/cmake/GlobalConfig.cmake)
include(../python/cmake/QifConfig.cmake)

add_subdirectory(pybind11)

# build qif python module
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/qif/*.cpp)
pybind11_add_module(qif_module qif.cpp ${SOURCES})

set_target_properties(qif_module PROPERTIES OUTPUT_NAME "qif")
target_link_libraries(qif_module PRIVATE qif)
target_include_directories(qif_module PRIVATE include)

# enable precompiled headers in cmake >= 3.16.0, EXCEPT on linux/clang
if(NOT SYSTEM_QIF AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0" AND NOT (UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
	target_precompile_headers(qif_module REUSE_FROM qif)
endif()

# run test script
add_custom_target(python_test DEPENDS python)
add_custom_command(TARGET python_test COMMAND PYTHONPATH=. python3 ../../python/tests/test.py)
