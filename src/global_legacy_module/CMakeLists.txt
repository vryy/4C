baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    ale
    art_net
    beam3
    beaminteraction
    bele
    binstrategy
    config
    constraint
    contact
    core
    elemag
    fluid
    fluid_ele
    fluid_xfluid
    lib
    lubrication_ele
    mat
    membrane
    module_registry
    mortar
    nurbs_discret
    particle_engine
    porofluidmultiphase_ele
    poromultiphase_scatra
    red_airways
    rigidsphere
    scatra_ele
    shell7p
    so3
    solid_ele
    solid_poro_ele
    structure_new
    thermo
    torsion3
    truss3
    w1
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
