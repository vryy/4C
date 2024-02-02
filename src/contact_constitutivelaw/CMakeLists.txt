baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    global_data
    inpar
    io
    lib
    mat
    )

if(BACI_WITH_MIRCO)
  list(APPEND _dependencies mirco::mirco_lib)
endif(BACI_WITH_MIRCO)

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
