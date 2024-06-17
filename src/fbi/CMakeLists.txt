four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    beam3
    beaminteraction
    config
    core
    fluid
    fluid_ele
    geometry_pair
    global_data
    inpar
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
