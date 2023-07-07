baci_auto_define_module()

set(_beaminteraction_dependencies
    # cmake-format: sortable
    beam3
    beamcontact
    binstrategy
    core
    geometry_pair
    inpar
    io
    lib
    linalg
    mat
    rebalance
    rigidsphere
    so3
    structure
    structure_new
    timestepping
    truss3
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_beaminteraction_dependencies})
