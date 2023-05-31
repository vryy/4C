baci_auto_define_module()

set(_cut_dependencies
    # cmake-format: sortable
    core
    geometry
    inpar
    io
    lib
    linalg
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_cut_dependencies})
