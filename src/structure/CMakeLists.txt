baci_auto_define_module()

set(_structure_dependencies
    # cmake-format: sortable
    adapter
    beamcontact
    beaminteraction
    cardiovascular0d
    comm
    constraint
    contact
    discsh3
    inpar
    io
    lib
    linalg
    mat
    mor
    mortar
    poroelast
    so3
    stru_multi
    structure_new
    surfstress
    timestepping
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_structure_dependencies})
