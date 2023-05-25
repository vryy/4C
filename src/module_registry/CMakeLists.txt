baci_auto_define_module()

set(_module_registry_dependencies
    # cmake-format: sortable
    ale
    art_net
    beam3
    beaminteraction
    bele
    binstrategy
    constraint
    contact
    discsh3
    elemag
    fluid_ele
    immersed_problem
    inv_analysis
    lubrication_ele
    mat
    membrane
    mortar
    nurbs_discret
    particle_engine
    porofluidmultiphase_ele
    red_airways
    rigidsphere
    s8
    scatra_ele
    so3
    thermo
    torsion3
    truss3
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_module_registry_dependencies})
