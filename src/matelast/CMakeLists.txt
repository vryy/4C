baci_auto_define_module()

set(_matelast_dependencies
    # cmake-format: sortable
    comm
    core
    inpar
    io
    lib
    mat
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_matelast_dependencies})
