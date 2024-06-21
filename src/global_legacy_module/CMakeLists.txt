four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    ale
    art_net
    beam3
    beaminteraction
    bele
    config
    constraint
    contact
    core
    elemag
    fluid
    fluid_ele
    fluid_xfluid
    lubrication_ele
    mat
    membrane
    module_registry
    mortar
    particle_engine
    porofluidmultiphase_ele
    poromultiphase_scatra
    red_airways
    rigidsphere
    scatra_ele
    shell7p
    so3
    solid_3D_ele
    solid_poro_3D_ele
    solid_scatra_3D_ele
    structure_new
    thermo
    torsion3
    truss3
    w1
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
