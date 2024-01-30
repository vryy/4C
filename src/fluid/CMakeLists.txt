baci_auto_define_module()

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
    io
    lib
    mat
    mortar
    nurbs_discret
    poroelast
    red_airways
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
