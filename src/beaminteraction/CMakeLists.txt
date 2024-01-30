baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    beam3
    beamcontact
    binstrategy
    config
    core
    geometry_pair
    global_data
    inpar
    io
    lib
    mat
    rigidsphere
    so3
    solid_ele
    structure
    structure_new
    timestepping
    truss3
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
