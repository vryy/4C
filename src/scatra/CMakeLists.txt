four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    binstrategy
    config
    core
    fluid
    fluid_turbulence
    global_data
    inpar
    io
    lib
    mat
    mortar
    nurbs_discret
    poromultiphase_scatra
    scatra_ele
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
