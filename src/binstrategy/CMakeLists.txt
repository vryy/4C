baci_auto_define_module()

set(_binstrategy_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    bele
    comm
    fluid_ele
    inpar
    io
    lib
    linalg
    mortar
    rebalance
    rigidsphere
    scatra_ele
    so3
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_binstrategy_dependencies})
