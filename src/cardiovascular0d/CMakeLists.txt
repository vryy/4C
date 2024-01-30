baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    global_data
    inpar
    io
    lib
    mor
    so3
    solver_nonlin_nox
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
