four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    binstrategy
    config
    contact
    core
    global_data
    inpar
    lib
    mat
    nurbs_discret
    poroelast
    poroelast_scatra
    so3
    solver_nonlin_nox
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
