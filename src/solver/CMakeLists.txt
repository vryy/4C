baci_auto_define_module()

set(_solver_dependencies
    # cmake-format: sortable
    inpar
    io
    lib
    linalg
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_solver_dependencies})
