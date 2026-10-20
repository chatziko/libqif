include_guard(GLOBAL)

# use backward.cpp of pretty stacktraces
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../../external/backward-cpp/CMakeLists.txt")
	set(CMAKE_UNITY_BUILD OFF)									# never use unity build for backward.cpp, it breaks
	add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../../external/backward-cpp backward-cpp)
	set(CMAKE_UNITY_BUILD ${DEFAULT_UNITY_BUILD})

else()
	# backward.cpp not present, add dummy commands to avoid errors
	set(BACKWARD_ENABLE "")
	macro(add_backward target)
	endmacro()
endif()
