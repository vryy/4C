baci_auto_define_module()

set(_so3_dependencies
    # cmake-format: sortable
    comm
    contact
    core
    fiber
    fluid_ele
    immersed_problem
    inpar
    io
    lib
    linalg
    mat
    nurbs_discret
    structure_new
    surfstress
    thermo
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_so3_dependencies})
