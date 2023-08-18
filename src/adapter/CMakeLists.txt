baci_auto_define_module()

set(_adapter_dependencies
    # cmake-format: sortable
    ale
    beaminteraction
    comm
    constraint
    contact
    core
    fluid
    fpsi
    fsi
    io
    levelset
    lib
    lubrication
    mat
    membrane
    mortar
    nurbs_discret
    pasi
    poroelast
    porofluidmultiphase
    s8
    scatra
    so3
    solid_ele
    ssi
    struct_ale
    structure
    structure_new
    thermo
    xfem
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_adapter_dependencies})
