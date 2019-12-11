cmake_minimum_required(VERSION 2.6)

project(libqif)

add_subdirectory(../external/pybind11 pybind11)

add_custom_target(python)	# target that build all python bindings

add_subdirectory(qif)

# run test script
add_custom_target(python_test DEPENDS python)
add_custom_command(TARGET python_test COMMAND PYTHONPATH=. python3 ../../python/tests/test.py)
