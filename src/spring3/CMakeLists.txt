baci_auto_define_module()

set(_spring3_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    core
    inpar
    lib
    mat
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_spring3_dependencies})
