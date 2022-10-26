#! This deprecated function links all libraries to a given target.#
#
function(deprecated_baci_add_library_linked_with_everything target)
  add_library(${target} ${ARGN})
  target_link_libraries(${target} PUBLIC ${BACI_ALL_ENABLED_EXTERNAL_LIBS})
endfunction()
