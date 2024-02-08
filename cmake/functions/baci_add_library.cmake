#
# This function adds a library target with the given sources and headers and sets all BACI-wide defaults on it.
# The type of target is determined automatically depending on whether any sources or headers are given,
#
# Usage:
#   baci_add_library(<target_name>
#     [SOURCES <list of source files>]
#     [HEADERS <list of header files>])
#
function(baci_add_library _target)
  # Parse arguments
  set(options "")
  set(oneValueArgs "")
  set(multiValueArgs SOURCES HEADERS)
  cmake_parse_arguments(
    _parsed
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if(DEFINED _parsed_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "There are unparsed arguments: ${_parsed_UNPARSED_ARGUMENTS}")
  endif()

  # Define an interface library for usage requirements only
  add_library(${_target}_deps INTERFACE)

  # See if there are any sources and headers and define appropriate targets.
  if(NOT _parsed_SOURCES)
    if(_parsed_HEADERS)
      message(DEBUG "${_target} is a header-only target")
    else()
      message(DEBUG "${_target} is a pure interface target without headers or sources")
    endif()
  else()
    message(DEBUG "${_target} is a target with sources")
    add_library(${_target}_objs OBJECT ${_parsed_SOURCES})
    target_link_libraries(${_target}_objs PUBLIC ${_target}_deps)
  endif()

  # Check that every header includes baci_config.H
  foreach(_header ${_parsed_HEADERS})
    file(READ ${_header} _header_content)
    string(FIND "${_header_content}" "#include \"baci_config.H\"" _index)
    if(_index EQUAL -1)
      message(
        FATAL_ERROR
          "The file \"${_header}\" does not include \"baci_config.H\". Please include the file and rerun CMake."
        )
    endif()
  endforeach()

  # Add the headers as a file set, which will automatically add the current directory as a search directory
  target_sources(${_target}_deps INTERFACE FILE_SET HEADERS FILES ${_parsed_HEADERS})

  # Add all global compile definitions
  target_link_libraries(${_target}_deps INTERFACE baci_global_compile_settings)

  # Link against all default external libraries
  baci_link_default_external_libraries(${_target}_deps INTERFACE)
endfunction()
