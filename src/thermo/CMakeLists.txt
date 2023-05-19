baci_auto_define_module()

set(_thermo_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    geometry
    inpar
    io
    lib
    linalg
    mat
    nurbs_discret
    solver
    timestepping
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_thermo_dependencies})
