include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.13)

include(../misc/cmake/GlobalHead.cmake)
project(libqif)
include(../misc/cmake/GlobalConfig.cmake)

# If this dir is the root, we need to include qif (if not qif is already loaded, and add_subdirectory will do nothing)
add_subdirectory(../lib lib)

# build pybind11
set(PYBIND11_PYTHON_VERSION 3.6)	# require at least python 3.6
add_subdirectory(../external/pybind11 pybind11)

# build the python module inside the "qif" dir, unless CMAKE_LIBRARY_OUTPUT_DIRECTORY is set (see setup.py)
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/qif")
endif()

# build qif python module
file(GLOB_RECURSE CPP_SOURCES src/*.cpp)
pybind11_add_module(qif_module ${CPP_SOURCES})

set_target_properties(qif_module PROPERTIES OUTPUT_NAME "_qif")		# the cpp module is build as "_qif", then imported by __init__.py
target_link_libraries(qif_module PRIVATE qif)
target_include_directories(qif_module PRIVATE include)

# enable precompiled headers in cmake >= 3.16.0, EXCEPT on linux/clang
if("${BUILD_QIF}" AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0" AND NOT (UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
	target_precompile_headers(qif_module REUSE_FROM qif)
endif()

# copy .py files to the output dir (preserving directories)
file(COPY src/qif DESTINATION . FILES_MATCHING PATTERN "*.py")

# on windows, copy required dlls (mp++.dll and libopenblas.dll) to the target dir
if(MSVC)
	find_file(OPENBLAS_PATH NAME libopenblas.dll PATH_SUFFIXES lib/ REQUIRED)
	find_file(MSVPC140_PATH NAME msvcp140.dll REQUIRED)

	add_custom_command(
		TARGET qif_module POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy
				$<TARGET_FILE:mp++> ${OPENBLAS_PATH} ${MSVPC140_PATH}
				$<TARGET_FILE_DIR:qif_module>)
endif()

# run test script
# add_custom_target(python_test DEPENDS python)
# add_custom_command(TARGET python_test COMMAND PYTHONPATH=. python3 ../../python/tests/test.py)
