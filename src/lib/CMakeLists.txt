baci_auto_define_module()

set(_lib_dependencies
    # cmake-format: sortable
    ale
    beaminteraction
    binstrategy
    comm
    contact_constitutivelaw
    core
    fiber
    fluid
    fluid_ele
    fluid_xfluid
    inpar
    io
    mat
    nurbs_discret
    particle_engine
    poromultiphase_scatra
    rebalance
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_lib_dependencies})
