baci_auto_define_module()

set(_mortar_dependencies
    # cmake-format: sortable
    binstrategy
    contact
    core
    inpar
    io
    lib
    mat
    nurbs_discret
    poroelast
    rebalance
    so3
    solver_nonlin_nox
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_mortar_dependencies})
