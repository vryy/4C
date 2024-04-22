# Declare an option with given name, description and default value (ON or OFF).
# This function is almost equivalent to the builtin CMake option() command,
# except that it also prints info on whether the option is set.
# Options need to start with "BACI_".
function(baci_process_global_option option_name description default)
  if(NOT option_name MATCHES "(BACI_|FOUR_C_).*")
    message(
      FATAL_ERROR
        "Disallowed option '${option_name}'. Option needs to start with 'BACI_' or 'FOUR_C_'."
      )
  endif()
  option("${option_name}" "${description}" "${default}")
  if(${option_name})
    message(STATUS "Option ${option_name} = ON")
  else()
    message(STATUS "Option ${option_name} = OFF")
  endif()
endfunction()
