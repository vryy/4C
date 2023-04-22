# This function adds an internal library to a module. The libary will be named <module>_<name>.
# Parameters:
#   MODULE name of the module which the library belongs to
#   NAME name of the library to add
#   SOURCES list of all source files
#   DEPENDENCIES list of all libraries that are required by the newly defined library
function(baci_add_library_to_module)
  set(options "")
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

  if(DEFINED parsed_KEYWORDS_MISSING_VALUES)
    message(SEND_ERROR "There are missing values for keywords: ${parsed_KEYWORDS_MISSING_VALUES}")
  endif()

  set(target_name "${parsed_MODULE}_${parsed_NAME}")
  baci_add_library(${target_name} ${parsed_SOURCES})
  baci_add_dependency(${target_name} ${parsed_DEPENDENCIES})
  target_link_libraries(${parsed_MODULE} INTERFACE ${target_name})
endfunction()
