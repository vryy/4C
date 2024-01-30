baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beaminteraction
    config
    core
    global_data
    inpar
    io
    lib
    mat
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
