four_c_auto_define_module()

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
    mat
    mortar
    poroelast
    so3
    solver_nonlin_nox
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
