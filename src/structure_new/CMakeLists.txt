baci_auto_define_module()

set(_structure_new_dependencies
    # cmake-format: sortable
    adapter
    beam3
    beamcontact
    beaminteraction
    browniandyn
    cardiovascular0d
    comm
    constraint
    contact
    contact_aug
    core
    geometry
    inpar
    io
    lib
    linalg
    loca_continuation
    mat
    so3
    solver
    solver_nonlin_nox
    struct_ale
    timestepping
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_structure_new_dependencies})
