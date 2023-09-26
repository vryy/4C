baci_auto_define_module()

set(_dependencies
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
    lubrication_ele
    mat
    membrane
    module_registry
    mortar
    nurbs_discret
    particle_engine
    porofluidmultiphase_ele
    red_airways
    rigidsphere
    scatra_ele
    shell7p
    so3
    thermo
    torsion3
    truss3
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
