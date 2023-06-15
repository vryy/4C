baci_auto_define_module()

set(_surfstress_dependencies
    # cmake-format: sortable
    inpar
    io
    lib
    linalg
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_surfstress_dependencies})
