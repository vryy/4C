baci_auto_define_module()

set(_bele_dependencies
    # cmake-format: sortable
    core
    lib
    mat
    so3
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_bele_dependencies})
