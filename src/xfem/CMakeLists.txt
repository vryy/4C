baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    bele
    config
    contact
    core
    cut
    fluid_ele
    global_data
    inpar
    io
    lib
    mat
    mortar
    poroelast
    so3
    solver_nonlin_nox
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
