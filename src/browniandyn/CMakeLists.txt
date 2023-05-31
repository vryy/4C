baci_auto_define_module()

set(_browniandyn_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    io
    lib
    linalg
    rigidsphere
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_browniandyn_dependencies})
