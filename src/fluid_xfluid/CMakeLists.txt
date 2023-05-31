baci_auto_define_module()

set(_fluid_xfluid_dependencies
    # cmake-format: sortable
    adapter
    cut
    fluid
    fluid_ele
    inpar
    io
    levelset
    lib
    linalg
    mat
    solver
    xfem
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fluid_xfluid_dependencies})
