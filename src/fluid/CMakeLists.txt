four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    ale
    art_net
    config
    core
    fluid_ele
    fluid_turbulence
    global_data
    inpar
    mat
    mortar
    poroelast
    red_airways
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
