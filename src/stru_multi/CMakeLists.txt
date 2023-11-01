baci_auto_define_module()

set(_stru_multi_dependencies
    # cmake-format: sortable
    inpar
    io
    lib
    mat
    so3
    structure
    surfstress
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_stru_multi_dependencies})
