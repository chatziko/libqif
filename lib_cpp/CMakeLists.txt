include_guard(GLOBAL)
cmake_minimum_required(VERSION 3.13)

project(libqif)

# create qif_cpp target that either builds shared lib, or uses the system one
#
if("${BUILD_QIF}")
    find_package(GSL REQUIRED)

	# build armadillo in dir 'armadillo'
	set(CMAKE_UNITY_BUILD OFF)								# never use unity for armadillo, it breaks
	add_subdirectory(../external/armadillo armadillo EXCLUDE_FROM_ALL)
	set(CMAKE_UNITY_BUILD ${DEFAULT_UNITY_BUILD})

	# build mp++ in dir 'mppp'
	add_subdirectory(../external/mppp mppp EXCLUDE_FROM_ALL)

	# build osqp in dir 'osqp', don't include targets in 'all' (we'll only depend on osqpstatic)
	set(CMAKE_UNITY_BUILD OFF)								# never use unity for osqp, it breaks
	add_subdirectory(../external/osqp osqp EXCLUDE_FROM_ALL)
	set(CMAKE_UNITY_BUILD ${DEFAULT_UNITY_BUILD})

	file(GLOB LIB_SOURCES src/*.cpp)								# get all src/*.cpp files in LIB_SOURCES
	add_library(qif_cpp ${LIB_SOURCES})								# libqif_cpp, depends on all src/*.cpp files. SHARED/STATIC is decided based on BUILD_SHARED_LIBS (see GlobalConfig.cmake)
	set_target_properties(qif_cpp PROPERTIES OUTPUT_NAME "qif")		# the target is named "qif_cpp", the library is named just "qif"

	target_include_directories(qif_cpp PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include" include)

else()
	# using existing qif
	find_library(QIF_PATH NAME qif PATH_SUFFIXES lib/ REQUIRED)

	add_library(qif_cpp INTERFACE)
	target_link_libraries(qif_cpp INTERFACE ${QIF_PATH})							# targets using qif_cpp will link to the system qif
endif()

# linked libraries
#
# PRIVATE, only needed for qif_cpp itself
if("${BUILD_QIF}")
	target_link_libraries(qif_cpp PRIVATE osqpstatic GSL::gsl armadillo mp++)	# armadillo/mp++ are not really needed, we only list then to get their include_directories

	if(${QIF_USE_ORTOOLS})
		target_link_libraries(qif_cpp PRIVATE ortools)							# also listed to get its include_directories
	endif()
	if(${QIF_USE_GLPK})
		# we use a wrapper so only qif_cpp itself needs to be linked with glpk
		target_link_libraries(qif_cpp PRIVATE ${LIB_GLPK})
	endif()
endif()

# INTERFACE, needed for targets _using_ qif_cpp
target_link_libraries(qif_cpp INTERFACE armadillo mp++)
if(${QIF_USE_ORTOOLS})
	target_link_libraries(qif_cpp INTERFACE ortools)
endif()
if(UNIX AND NOT APPLE)
	# dl and rt are needed in manylinux
	target_link_libraries(qif_cpp INTERFACE ${CMAKE_DL_LIBS} rt)
endif()

# compiler options
# * we build for march=sandybridge by default, can be changed with -DMARCH=...
#
target_compile_options(qif_cpp INTERFACE -Wall)
if(MSVC)
	target_compile_options(qif_cpp INTERFACE /permissive- /arch:$<IF:$<BOOL:${ARCH}>,${ARCH},AVX>)
else()
	target_compile_options(qif_cpp INTERFACE -Wextra -pedantic -pthread -march=$<IF:$<BOOL:${MARCH}>,${MARCH},sandybridge>)
endif()
target_compile_definitions(qif_cpp INTERFACE $<$<CONFIG:Release>:ARMA_NO_DEBUG>)

if("${BUILD_QIF}")
	# if we're building qif_cpp, copy to PRIVATE so that they are used also for qif_cpp itself
	target_compile_options    (qif_cpp PRIVATE $<TARGET_PROPERTY:qif_cpp,INTERFACE_COMPILE_OPTIONS>)
	target_compile_definitions(qif_cpp PRIVATE $<TARGET_PROPERTY:qif_cpp,INTERFACE_COMPILE_DEFINITIONS>)
endif()

# config file
configure_file(include/qif_bits/config.h include/qif_bits/config.h)



# installation
#
if("${BUILD_QIF}")	# can only install if we built
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION "include"
		PATTERN "precompiled.h" EXCLUDE
		PATTERN "osqp" EXCLUDE
		PATTERN "miniball" EXCLUDE
		PATTERN "_*" EXCLUDE)
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION "include")

	# install armadillo's includes
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/armadillo/tmp/include/ DESTINATION "include")

	# install mp++'s includes, they have to be globally installed cause mp++ uses #include<...>
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../external/mppp/include/ DESTINATION "include")
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mppp/include/ DESTINATION "include")

	# install OSQP's includes under qif_bits/osqp
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../external/osqp/include/ DESTINATION "include/qif_bits/osqp")

	# install miniball's includes under qif_bits/miniball
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../external/miniball/cpp/main/ DESTINATION "include/qif_bits/miniball")

	install(TARGETS qif_cpp armadillo mp++ LIBRARY)
endif()

