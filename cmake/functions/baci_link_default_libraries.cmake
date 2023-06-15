# This function links all default external dependencies to the given target.
# The keyword defines the usual linking visibility (PUBLIC, INTERFACE, PRIVATE)
function(baci_link_default_external_libraries target keyword)
  target_link_libraries(${target} ${keyword} ${BACI_ALL_ENABLED_EXTERNAL_LIBS})
endfunction()
