# Automatically create a library for the sources and headers in the current directory. The target
# will be named based on the folder name and any parent modules. This name is returned in the variable
# AUTO_DEFINED_MODULE_NAME which is set at the call site.
function(baci_auto_define_module)
  # Get a name based on the current folder
  get_filename_component(_target ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  # Check if we are currently defining a part of a parent module. In this case, concatenate the module names
  if(NOT "${BACI_CURRENTLY_DEFINED_PARENT_MODULE}" STREQUAL "")
    set(_target "${BACI_CURRENTLY_DEFINED_PARENT_MODULE}_${_target}")
  endif()

  message(VERBOSE "Defining BACI module target: ${_target}")

  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  file(
    GLOB _sources
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.cpp *.cc *.c
    )
  file(
    GLOB _headers
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.h *.hpp
    )
  # Remove headers that only contain template instantiations
  list(
    FILTER
    _headers
    EXCLUDE
    REGEX
    "_fwd\.h(pp)?|\.inst\.[hH]"
    )

  baci_add_library(
    ${_target}
    SOURCES
    ${_sources}
    HEADERS
    ${_headers}
    )

  # Check if we are currently defining a part of a parent module and add ourselves as a dependency.
  if(NOT "${BACI_CURRENTLY_DEFINED_PARENT_MODULE}" STREQUAL "")
    baci_add_dependency(${BACI_CURRENTLY_DEFINED_PARENT_MODULE} ${_target})
  endif()

  # Now check if there are more directories that contain CMakeLists.txt. If yes, we also add those.
  # For this action, we become the parent module of the submodules we are about to define.
  set(BACI_CURRENTLY_DEFINED_PARENT_MODULE ${_target})

  # N.B. We need to directly glob for CMakeLists.txt files here to ensure
  # the glob reruns when a new CMakeLists.txt is added.
  file(
    GLOB children
    RELATIVE ${CMAKE_CURRENT_LIST_DIR}
    CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/*/CMakeLists.txt
    )
  foreach(child ${children})
    get_filename_component(_subdir ${child} DIRECTORY)
    add_subdirectory(${_subdir})
  endforeach()

  list(POP_BACK CMAKE_MESSAGE_INDENT)

  # Collect all modules in the global library
  target_link_libraries(${BACI_LIBRARY_NAME} PUBLIC ${_target}_deps)
  target_link_libraries(${BACI_LIBRARY_NAME} PRIVATE $<TARGET_NAME_IF_EXISTS:${_target}_objs>)

  # Simulate a "return" by setting a variable at the call site
  set(AUTO_DEFINED_MODULE_NAME
      ${_target}
      PARENT_SCOPE
      )
endfunction()
