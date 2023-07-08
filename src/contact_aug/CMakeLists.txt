baci_auto_define_module()

set(_contact_aug_dependencies
    # cmake-format: sortable
    contact
    core
    inpar
    io
    lib
    mortar
    solver_nonlin_nox
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_contact_aug_dependencies})
