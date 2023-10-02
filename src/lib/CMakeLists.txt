baci_auto_define_module()

set(_lib_dependencies
    # cmake-format: sortable
    beaminteraction
    binstrategy
    comm
    contact_constitutivelaw
    core
    fiber
    fluid
    fluid_xfluid
    global_legacy_module
    inpar
    io
    mat
    nurbs_discret
    particle_engine
    poromultiphase_scatra
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_lib_dependencies})
