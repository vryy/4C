baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    binstrategy
    config
    contact_aug
    contact_constitutivelaw
    core
    global_data
    inpar
    io
    lib
    mat
    mortar
    scatra
    scatra_ele
    so3
    solver_nonlin_nox
    structure_new
    xfem
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
