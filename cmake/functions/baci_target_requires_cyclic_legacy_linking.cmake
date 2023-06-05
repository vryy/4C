# Call this function on a target that needs to link to all modules
# This is a legacy call and should NOT be applied to new targets
function(baci_target_requires_cyclic_legacy_linking target)
  set_property(GLOBAL APPEND PROPERTY BACI_TARGETS_REQUIRING_CYCLIC_LEGACY_LINKING ${target})
endfunction()
