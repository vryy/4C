baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    io
    lib
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

option(
  DEBUGCUTLIBRARY "Turn on special output inside the cut library to be able to track errors" OFF
  )
if(DEBUGCUTLIBRARY)
  target_compile_definitions(${AUTO_DEFINED_MODULE_NAME} PUBLIC DEBUGCUTLIBRARY)
endif()
