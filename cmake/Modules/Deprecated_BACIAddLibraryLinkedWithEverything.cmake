#! This deprecated function creates a library and links all required dependencies to it.
#
function(deprecated_baci_add_library_linked_with_everything target)
  add_library(${target} ${ARGN})
  baci_link_default_libraries(${target})
endfunction()
