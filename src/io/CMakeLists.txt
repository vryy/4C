baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beam3
    config
    core
    fluid_ele
    global_data
    inpar
    lib
    mat
    nurbs_discret
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
