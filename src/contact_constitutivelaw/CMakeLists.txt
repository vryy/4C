four_c_auto_define_module()

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

if(FOUR_C_WITH_MIRCO)
  four_c_add_external_dependency(${AUTO_DEFINED_MODULE_NAME} mirco::mirco_lib)
endif(FOUR_C_WITH_MIRCO)

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
