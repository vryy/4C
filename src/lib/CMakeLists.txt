baci_auto_define_module()

set(_lib_dependencies
    # cmake-format: sortable
    ale
    beaminteraction
    binstrategy
    comm
    contact_constitutivelaw
    core
    create_rtdfiles
    elemag
    fiber
    fluid
    fluid_ele
    fluid_xfluid
    immersed_problem
    inpar
    io
    linalg
    mat
    nurbs_discret
    particle_engine
    poromultiphase_scatra
    rebalance
    scatra_ele
    structure_new
    xfem
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_lib_dependencies})
