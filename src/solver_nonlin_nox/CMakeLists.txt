baci_auto_define_module()

set(_solver_nonlin_nox_dependencies
    # cmake-format: sortable
    cardiovascular0d
    constraint
    contact
    contact_aug
    core
    inpar
    lib
    scatra
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_solver_nonlin_nox_dependencies})
