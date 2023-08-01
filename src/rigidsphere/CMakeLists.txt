baci_auto_define_module()

set(_rigidsphere_dependencies
    # cmake-format: sortable
    beaminteraction
    core
    inpar
    lib
    mat
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_rigidsphere_dependencies})
