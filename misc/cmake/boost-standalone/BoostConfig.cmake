# Minimal "Boost" package, used only by SoPlex (external/soplex), which calls find_package(Boost) but only needs
# Boost.Multiprecision. We provide it in standalone mode, using the external/boost-multiprecision and
# external/boost-config submodules, so that no system Boost is needed.

set(Boost_VERSION "1.92.0")
set(Boost_VERSION_STRING "1.92.0")
set(Boost_VERSION_MACRO 109200)
set(Boost_FOUND TRUE)

if(NOT TARGET Boost::boost)
	set(_qif_external "${CMAKE_CURRENT_LIST_DIR}/../../../external")

	add_library(Boost::boost INTERFACE IMPORTED GLOBAL)
	set_target_properties(Boost::boost PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${_qif_external}/boost-multiprecision/include;${_qif_external}/boost-config/include;${CMAKE_CURRENT_LIST_DIR}/include"
		INTERFACE_COMPILE_DEFINITIONS "BOOST_MP_STANDALONE")
	add_library(Boost::headers ALIAS Boost::boost)
endif()
