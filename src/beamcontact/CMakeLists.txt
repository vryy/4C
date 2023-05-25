baci_auto_define_module()

set(_beamcontact_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    contact
    core
    inpar
    io
    lib
    linalg
    rigidsphere
    structure
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_beamcontact_dependencies})
