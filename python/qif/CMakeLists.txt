cmake_minimum_required(VERSION 2.6)

project(libqif)

pybind11_add_module(channel channel.cpp)

target_link_libraries(channel LINK_PRIVATE qif)

add_dependencies(python channel)