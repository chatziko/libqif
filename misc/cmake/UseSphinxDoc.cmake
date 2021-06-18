find_package(PythonInterp REQUIRED)
find_package(Sphinx REQUIRED)

function(add_sphinx_doc)
  set(options)
  set(oneValueArgs
    SOURCE_DIR
    PYTHONPATH
    TARGET_NAME
  )
  set(multiValueArgs)

  cmake_parse_arguments(SPHINX_DOC
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )

  add_custom_target(${SPHINX_DOC_TARGET_NAME}
    # copy sources to _build/sources so that autosummary files are generated there and not in the source dir
    COMMAND
      rm -rf ${CMAKE_CURRENT_BINARY_DIR}/_build/ &&
      mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/_build/ &&
      cp -r ${SPHINX_DOC_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/_build/sources &&
      PYTHONPATH=${SPHINX_DOC_PYTHONPATH} ${SPHINX_EXECUTABLE}
         -q
         -b html
         -d ${CMAKE_CURRENT_BINARY_DIR}/_build/doctrees
         ${CMAKE_CURRENT_BINARY_DIR}/_build/sources
         ${CMAKE_CURRENT_BINARY_DIR}/_build/html
    VERBATIM
  )

  message(STATUS "Added ${SPHINX_DOC_TARGET_NAME} [Sphinx] target to build documentation")
endfunction()
