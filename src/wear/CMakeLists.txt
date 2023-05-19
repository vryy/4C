baci_auto_define_module()

set(_wear_dependencies
    # cmake-format: sortable
    adapter
    ale
    contact
    contact_aug
    core
    fs3i
    fsi
    inpar
    lib
    linalg
    mortar
    nurbs_discret
    so3
    solver
    structure
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_wear_dependencies})
