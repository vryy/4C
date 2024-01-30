baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    io
    solver_nonlin_nox
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
