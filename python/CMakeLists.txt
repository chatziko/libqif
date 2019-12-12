cmake_minimum_required(VERSION 2.6)

project(libqif)

include(cmake/Common.cmake)

# if the bindings are compiled alone (not within a library build), add an interface qif target
# that has the same options as the normal qif target and links to the system-wide qif.
if(NOT TARGET qif)
	find_library(SYSTEM_QIF NAME qif PATH_SUFFIXES lib/)

	add_library(qif INTERFACE)
	target_link_libraries(qif INTERFACE ${SYSTEM_QIF} ${QIF_LIBS})	# libraries to link

	target_compile_options(qif INTERFACE -Wall -Wextra -pedantic  $<$<CONFIG:Release>:-march=native> -fPIC -pthread)
	target_compile_definitions(qif INTERFACE $<$<CONFIG:Release>:ARMA_NO_DEBUG> qif_EXPORTS)
endif()

add_subdirectory(pybind11)

include_directories(include)

add_custom_target(python)	# target that build all python bindings

add_subdirectory(qif)

# run test script
add_custom_target(python_test DEPENDS python)
add_custom_command(TARGET python_test COMMAND PYTHONPATH=. python3 ../../python/tests/test.py)
