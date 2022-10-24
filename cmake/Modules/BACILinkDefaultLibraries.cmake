# This function links all default dependencies to the given target.
# This includes all external libraries and the drt_lib_error target.
function(baci_link_default_libraries target)
  target_link_libraries(${target} PUBLIC ${BACI_ALL_ENABLED_EXTERNAL_LIBS})
  target_link_libraries(${target} PUBLIC drt_lib_error)
endfunction()
