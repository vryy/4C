baci_auto_define_module()

set(_porofluidmultiphase_dependencies
    # cmake-format: sortable
    adapter
    art_net
    inpar
    io
    lib
    mat
    porofluidmultiphase_ele
    poromultiphase_scatra
    rebalance
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_porofluidmultiphase_dependencies})
