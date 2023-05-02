# This function links all default dependencies to the given target.
# This includes all external libraries and the drt_lib_error target.
# The keyword defines the usual linking visibility (PUBLIC, INTERFACE, PRIVATE)
function(baci_link_default_libraries target keyword)
  target_link_libraries(${target} ${keyword} ${BACI_ALL_ENABLED_EXTERNAL_LIBS})
  target_link_libraries(${target} ${keyword} lib_error)
endfunction()
