four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    config
    core
    global_data
    inpar
    io
    lib
    so3
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
