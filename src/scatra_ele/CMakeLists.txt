baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    fiber
    fluid
    fluid_ele
    global_data
    inpar
    io
    lib
    mat
    nurbs_discret
    porofluidmultiphase_ele
    so3
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
