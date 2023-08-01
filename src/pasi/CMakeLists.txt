baci_auto_define_module()

set(_pasi_dependencies
    # cmake-format: sortable
    adapter
    comm
    inpar
    io
    lib
    particle_algorithm
    particle_wall
    solver_nonlin_nox
    structure
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_pasi_dependencies})
