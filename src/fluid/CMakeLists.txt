baci_auto_define_module()

set(_fluid_dependencies
    # cmake-format: sortable
    adapter
    art_net
    comm
    core
    fluid_ele
    fluid_turbulence
    inpar
    io
    levelset
    lib
    linalg
    mat
    mortar
    nurbs_discret
    poroelast
    rebalance
    red_airways
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fluid_dependencies})
