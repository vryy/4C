function(baci_link_executable_with_all executable)
  if(NOT BACI_ALL_DEFINED_MODULE_TARGETS)
    message(SEND_ERROR "Tried to link an executable before all modules were defined.")
  endif()

  target_link_libraries(
    ${executable} PRIVATE -Wl,--start-group ${BACI_ALL_DEFINED_MODULE_TARGETS} -Wl,--end-group
    )
endfunction()
