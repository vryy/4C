# Declare an option with given name, description and default value (ON or OFF).
# This function is almost equivalent to the builtin CMake option() command,
# except that it also checks whether the option is set. In this case, the option is
# appended to the global property BACI_GLOBAL_COMPILE_DEFINITIONS. For later
# consumption by baci_add_library().
function(baci_process_global_option option_name description default)
  option("${option_name}" "${description}" "${default}")
  if(${option_name})
    message(STATUS "Option ${option_name} = ON")
    set_property(GLOBAL APPEND PROPERTY BACI_GLOBAL_COMPILE_DEFINITIONS "${option_name}")
  else()
    message(STATUS "Option ${option_name} = OFF")
  endif()
endfunction()
