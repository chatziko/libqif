include_guard(GLOBAL)

# Global config, included in all (potentially) top-level CMakeLists: qif / samples / python
# before the project() declaration

# Allow setting some options from the environment
foreach (VAR CMAKE_TOOLCHAIN_FILE CMAKE_INSTALL_PREFIX BUILD_QIF)
	if(DEFINED ENV{${VAR}})
		set(${VAR} "$ENV{${VAR}}" CACHE PATH "" FORCE)
	endif()
endforeach()