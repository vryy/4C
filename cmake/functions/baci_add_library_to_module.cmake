# This function adds an internal library to a module. The libary will be named <module>_<name>.
# Parameters:
#   MODULE name of the module which the library belongs to
#   NAME name of the library to add
#   SOURCES list of all source files
#   DEPENDENCIES list of all libraries that are required by the newly defined library
function(baci_add_library_to_module)
  set(options INTERFACE)
  set(oneValueArgs MODULE NAME)
  set(multiValueArgs SOURCES DEPENDENCIES)
  cmake_parse_arguments(
    parsed
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if(DEFINED parsed_UNPARSED_ARGUMENTS)
    message(SEND_ERROR "There are unparsed arguments: ${parsed_UNPARSED_ARGUMENTS}")
  endif()

  set(target_name "${parsed_MODULE}_${parsed_NAME}")

  if(parsed_INTERFACE)
    add_library(${target_name} INTERFACE)
    target_include_directories(
      ${target_name} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      )
    # link with external libraries
    baci_link_default_libraries(${target_name} INTERFACE)
  else()
    baci_add_library(${target_name} ${parsed_SOURCES})
  endif()

  if(parsed_DEPENDENCIES)
    baci_add_dependency(${target_name} ${parsed_DEPENDENCIES})
  endif()

  target_link_libraries(${parsed_MODULE} INTERFACE ${target_name})
endfunction()
