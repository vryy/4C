#! This deprecated function creates a library and links all required dependencies to it.
#
function(deprecated_baci_add_library_linked_with_everything target)
  add_library(${target} ${ARGN})

  # create a local copy of all targets and remove the target we currently define since target_link_libraries errors
  set(_baci_libraries_local_copy ${BACI_LIBRARIES})
  list(REMOVE_ITEM _baci_libraries_local_copy ${target})
  # Allow usage of all internal libraries
  target_link_libraries(${target} PUBLIC ${_baci_libraries_local_copy})

  # link external libraries
  baci_link_default_libraries(${target})

  # For legacy reasons and due to cyclic dependencies, we need to repeat libraries on the link line. The number of times
  # specified here is the minimum number that was necessary to build all executable targets in this project.
  set_target_properties(${target} PROPERTIES LINK_INTERFACE_MULTIPLICITY 7)
endfunction()
