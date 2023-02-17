function(baci_add_dependency target)
  target_link_libraries(${target} PUBLIC ${ARGN})
endfunction()
