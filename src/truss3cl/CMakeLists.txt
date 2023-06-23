baci_auto_define_module()

set(_truss3cl_dependencies
    # cmake-format: sortable
    core
    inpar
    lib
    linalg
    mat
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_truss3cl_dependencies})
