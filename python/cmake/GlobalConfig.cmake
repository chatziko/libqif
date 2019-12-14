include_guard(GLOBAL)

# Global config, included in all (potentially) top-level CMakeLists: qif / samples / python

# default type is Release, change with cmake -DCMAKE_BUILD_TYPE=
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# use unity builds by default for all targets
#
if(NOT CMAKE_UNITY_BUILD AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0") 
	set(CMAKE_UNITY_BUILD ON)
	set(DEFAULT_UNITY_BUILD ON)
	# set(CMAKE_UNITY_BUILD_BATCH_SIZE 0)
endif()

# use c++17
set(CMAKE_CXX_STANDARD 17)												# easier method for setting c++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
if(APPLE)
	LINK_DIRECTORIES(/usr/local/lib)
endif()

# Macros
MACRO(SUBDIRLIST result curdir)
	FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
	SET(dirlist "")
	FOREACH(child ${children})
		IF(IS_DIRECTORY ${curdir}/${child})
			LIST(APPEND dirlist ${child})
		ENDIF()
	ENDFOREACH()
	SET(${result} ${dirlist})
ENDMACRO()