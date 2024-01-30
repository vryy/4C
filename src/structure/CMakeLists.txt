baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    beamcontact
    beaminteraction
    cardiovascular0d
    config
    constraint
    contact
    core
    global_data
    inpar
    io
    lib
    mat
    mor
    mortar
    poroelast
    so3
    stru_multi
    structure_new
    timestepping
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
