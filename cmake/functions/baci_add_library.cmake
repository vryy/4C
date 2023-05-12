# This function adds a library target with the given sources.
function(baci_add_library target)
  add_library(${target} ${ARGN})
  # only allow header inclusions from the current source folder
  target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

  # link with external libraries
  baci_link_default_libraries(${target} PUBLIC)
endfunction()
