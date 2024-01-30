baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    global_data
    inpar
    io
    lib
    mortar
    solver_nonlin_nox
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
