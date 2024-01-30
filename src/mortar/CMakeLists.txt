baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    binstrategy
    config
    contact
    core
    global_data
    inpar
    io
    lib
    mat
    nurbs_discret
    poroelast
    poroelast_scatra
    so3
    solver_nonlin_nox
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
