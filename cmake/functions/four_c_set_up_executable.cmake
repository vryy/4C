# Call this function on an executable of this project.
# The executable will link to the main library and use the internal compiler settings.
function(four_c_set_up_executable target)
  target_link_libraries(${target} PRIVATE ${FOUR_C_LIBRARY_NAME})
  target_link_libraries(${target} PRIVATE four_c_private_compile_interface)
endfunction()
