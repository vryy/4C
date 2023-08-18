baci_auto_define_module()

set(_truss3_dependencies
    # cmake-format: sortable
    core
    inpar
    lib
    mat
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_truss3_dependencies})
