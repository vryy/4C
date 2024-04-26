four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    inpar
    io
    lib
    mortar
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
