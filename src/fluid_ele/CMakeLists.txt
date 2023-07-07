baci_auto_define_module()

set(_fluid_ele_dependencies
    # cmake-format: sortable
    bele
    core
    cut
    fluid
    immersed_problem
    inpar
    io
    lib
    linalg
    mat
    nurbs_discret
    poroelast
    xfem
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fluid_ele_dependencies})
