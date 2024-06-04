four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    fluid
    fluid_ele
    global_data
    inpar
    lib
    mat
    nurbs_discret
    porofluidmultiphase_ele
    so3
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
