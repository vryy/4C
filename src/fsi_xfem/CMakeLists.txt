baci_auto_define_module()

set(_fsi_xfem_dependencies
    # cmake-format: sortable
    adapter
    constraint
    contact
    core
    fluid_xfluid
    fsi
    inpar
    io
    lib
    linalg
    poroelast
    solver
    structure
    xfem
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fsi_xfem_dependencies})
