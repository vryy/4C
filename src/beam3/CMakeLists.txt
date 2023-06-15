baci_auto_define_module()

set(_beam3_dependencies
    # cmake-format: sortable
    beaminteraction
    core
    inpar
    io
    lib
    linalg
    mat
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_beam3_dependencies})
