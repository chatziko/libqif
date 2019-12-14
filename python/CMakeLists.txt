cmake_minimum_required(VERSION 2.6)

project(libqif)

include(../python/cmake/GlobalConfig.cmake)
include(../python/cmake/QifConfig.cmake)

add_subdirectory(pybind11)

include_directories(include)

add_custom_target(python)	# target that build all python bindings

add_subdirectory(qif)

# run test script
add_custom_target(python_test DEPENDS python)
add_custom_command(TARGET python_test COMMAND PYTHONPATH=. python3 ../../python/tests/test.py)
