baci_auto_define_module()

set(_contact_dependencies
    # cmake-format: sortable
    binstrategy
    contact_aug
    contact_constitutivelaw
    core
    inpar
    io
    lib
    mat
    mortar
    rebalance
    scatra
    scatra_ele
    so3
    solver_nonlin_nox
    structure_new
    xfem
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_contact_dependencies})
