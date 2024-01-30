baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    art_net
    config
    core
    global_data
    inpar
    io
    lib
    mat
    porofluidmultiphase_ele
    poromultiphase_scatra
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
