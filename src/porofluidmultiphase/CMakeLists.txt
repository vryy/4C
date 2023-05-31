baci_auto_define_module()

set(_porofluidmultiphase_dependencies
    # cmake-format: sortable
    adapter
    art_net
    geometry
    inpar
    io
    lib
    linalg
    mat
    porofluidmultiphase_ele
    poromultiphase_scatra
    rebalance
    solver
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_porofluidmultiphase_dependencies})
