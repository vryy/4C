baci_auto_define_module()

set(_s8_dependencies
    # cmake-format: sortable
    inpar
    lib
    linalg
    mat
    matelast
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_s8_dependencies})
