function(_baci_internal_link_with_debug_message target link_type deps)
  message(DEBUG "Linking: ${target} <- ${deps} (${link_type})")
  target_link_libraries(${target} ${link_type} ${deps})
endfunction()

function(baci_add_dependency target)
  # Internal BACI target in the library
  if(TARGET ${target}_deps)
    foreach(_dep ${ARGN})
      _baci_internal_link_with_debug_message(${target}_deps INTERFACE ${_dep}_deps)
    endforeach()
    # Some other target
  else()
    get_target_property(_target_type ${target} TYPE)
    if(${_target_type} STREQUAL INTERFACE_LIBRARY)
      foreach(_dep ${ARGN})
        _baci_internal_link_with_debug_message(${target} INTERFACE ${_dep}_deps)
      endforeach()
    elseif(${_target_type} STREQUAL EXECUTABLE)
      # Currently the code base is not able to link without cycles. When linking an executable to a libary, the
      # transitive dependencies must not contain any cycles. For static libraries this can be ignored when the
      # libraries are repeated on the link line. However, for shared libraries there is no such workaround.
      # For this reason, we issue an error here (for now).
      message(
        FATAL_ERROR "Trying to link parts of the library to ${target} which is an exectuable. "
                    "Please link against the whole library."
        )
    else()
      message(FATAL_ERROR "Cannot add dependency to ${target} of type ${_target_type}.")
    endif()
  endif()
endfunction()

function(baci_add_external_dependency target)
  # Internal BACI target in the library: adjust name
  if(TARGET ${target}_deps)
    set(target ${target}_deps)
  endif()
  get_target_property(_target_type ${target} TYPE)
  if(${_target_type} STREQUAL INTERFACE_LIBRARY)
    target_link_libraries(${target} INTERFACE ${ARGN})
  else()
    target_link_libraries(${target} PUBLIC ${ARGN})
  endif()
endfunction()
