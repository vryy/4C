baci_auto_define_module()

set(_cut_dependencies
    # cmake-format: sortable
    core
    inpar
    io
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_cut_dependencies})

option(
  DEBUGCUTLIBRARY "Turn on special output inside the cut library to be able to track errors" OFF
  )
if(DEBUGCUTLIBRARY)
  target_compile_definitions(${AUTO_DEFINED_MODULE_NAME} PUBLIC DEBUGCUTLIBRARY)
endif()
