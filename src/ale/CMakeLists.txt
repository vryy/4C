baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    contact
    core
    fluid_ele
    global_data
    inpar
    io
    lib
    mat
    mortar
    nurbs_discret
    so3
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
