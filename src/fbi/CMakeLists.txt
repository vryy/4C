baci_auto_define_module()

set(_fbi_dependencies
    # cmake-format: sortable
    adapter
    beam3
    beaminteraction
    binstrategy
    core
    fluid
    fluid_ele
    geometry_pair
    inpar
    io
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fbi_dependencies})
