four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    ale
    config
    contact
    contact_aug
    core
    fs3i
    fsi
    global_data
    inpar
    lib
    mortar
    nurbs_discret
    so3
    structure
    w1
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
