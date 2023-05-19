baci_auto_define_module()

set(_inv_analysis_dependencies
    # cmake-format: sortable
    adapter
    comm
    inpar
    io
    lib
    linalg
    mat
    solver
    timestepping
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_inv_analysis_dependencies})
