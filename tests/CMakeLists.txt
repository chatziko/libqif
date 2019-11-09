cmake_minimum_required(VERSION 2.6)

project(libqif)


add_subdirectory(../external/googletest googletest)					# build googletest in dir 'googletest'

file(GLOB_RECURSE TEST_SOURCES *.cpp)								# get all *.cpp files in TEST_SOURCES

add_executable(run ${TEST_SOURCES})
target_link_libraries(run qif gtest)								# crete 'run' executable, link against libqif and googletest

# enable precompiled headers in cmake >= 3.16.0
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
	target_precompile_headers(run REUSE_FROM qif)
endif()


# for 'make tests'
add_custom_target(tests)
add_dependencies(tests run)
