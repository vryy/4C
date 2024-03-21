function(_get_target_property_or_empty _result_var _target _property)
  get_target_property(${_result_var} ${_target} ${_property})
  if(${_result_var})
    return(PROPAGATE ${_result_var})
  else()
    set(${_result_var} "")
    return(PROPAGATE ${_result_var})
  endif()
endfunction()

#
# A wrapper around CMake's try_compile including all compiler setup done up to this
# point. All settings that are checked are added at once, and only if they all work,
# they will be used.
#
function(baci_add_settings_if_compiles _result_var)
  # Parse arguments
  set(options "")
  set(oneValueArgs "")
  set(multiValueArgs COMPILE_OPTIONS COMPILE_DEFINITIONS LINK_OPTIONS LINK_LIBRARIES)
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

  set(_test_source_code
      "#include <iostream>\nint main() { std::cout << \"Hello, world!\" << std::endl; return 0; }"
      )

  _get_target_property_or_empty(
    _compile_definitions baci_private_compile_interface INTERFACE_COMPILE_DEFINITIONS
    )
  _get_target_property_or_empty(
    _compile_options baci_private_compile_interface INTERFACE_COMPILE_OPTIONS
    )
  _get_target_property_or_empty(
    _link_libraries baci_private_compile_interface INTERFACE_LINK_LIBRARIES
    )
  _get_target_property_or_empty(_link_options baci_private_compile_interface INTERFACE_LINK_OPTIONS)

  # Note: this one must be a space-separated string!
  string(JOIN " " CMAKE_REQUIRED_FLAGS ${_compile_options})
  string(
    JOIN
    " "
    CMAKE_REQUIRED_FLAGS
    "${CMAKE_REQUIRED_FLAGS}"
    ${_parsed_COMPILE_OPTIONS}
    )

  # These ones are ;-lists.
  set(CMAKE_REQUIRED_DEFINITIONS ${_compile_definitions} ${_parsed_COMPILE_DEFINITIONS})
  set(CMAKE_REQUIRED_LINK_OPTIONS ${_link_options} ${_parsed_LINK_OPTIONS})
  set(CMAKE_REQUIRED_LIBRARIES ${_link_libraries} ${_parsed_LINK_LIBRARIES})

  check_cxx_source_compiles("${_test_source_code}" ${_result_var})

  if(${_result_var})
    target_compile_options(baci_private_compile_interface INTERFACE ${_parsed_COMPILE_OPTIONS})
    target_compile_definitions(
      baci_private_compile_interface INTERFACE ${_parsed_COMPILE_DEFINITIONS}
      )
    target_link_options(baci_private_compile_interface INTERFACE ${_parsed_LINK_OPTIONS})
    target_link_libraries(baci_private_compile_interface INTERFACE ${_parsed_LINK_LIBRARIES})
  endif()

endfunction()
