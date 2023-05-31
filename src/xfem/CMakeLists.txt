baci_auto_define_module()

set(_xfem_dependencies
    # cmake-format: sortable
    bele
    contact
    core
    cut
    fluid_ele
    geometry
    inpar
    io
    lib
    linalg
    mat
    mortar
    poroelast
    rebalance
    so3
    solver_nonlin_nox
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_xfem_dependencies})
