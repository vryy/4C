baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    config
    core
    global_data
    io
    lib
    rigidsphere
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
