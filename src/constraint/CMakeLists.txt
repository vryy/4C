baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    contact
    core
    global_data
    inpar
    io
    lib
    solver_nonlin_nox
    structure_new
    truss3
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
