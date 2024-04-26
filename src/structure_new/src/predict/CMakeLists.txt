four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    inpar
    io
    solver_nonlin_nox
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
