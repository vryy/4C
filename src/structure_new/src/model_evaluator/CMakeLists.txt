four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beam3
    beamcontact
    beaminteraction
    browniandyn
    cardiovascular0d
    config
    constraint
    constraint_framework
    contact
    contact_aug
    core
    global_data
    inpar
    solver_nonlin_nox
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
