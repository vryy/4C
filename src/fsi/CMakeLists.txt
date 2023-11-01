baci_auto_define_module()

set(_fsi_dependencies
    # cmake-format: sortable
    adapter
    ale
    binstrategy
    constraint
    core
    fbi
    fluid
    fluid_xfluid
    fsi_xfem
    inpar
    io
    lib
    mat
    mortar
    poroelast
    so3
    solver_nonlin_nox
    structure
    structure_new
    timestepping
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fsi_dependencies})
