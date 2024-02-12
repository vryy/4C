# Call this function on a target that needs to link to all modules
# This is a legacy call and should NOT be applied to new targets
function(baci_target_requires_cyclic_legacy_linking target)
  target_link_libraries(${target} PRIVATE ${BACI_LIBRARY_NAME})
endfunction()
