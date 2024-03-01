baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beam3
    beaminteraction
    bele
    config
    core
    fluid_ele
    global_data
    inpar
    io
    lib
    mortar
    rigidsphere
    scatra_ele
    so3
    solid_3D_ele
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
