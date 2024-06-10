four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    art_net
    config
    core
    global_data
    inpar
    mat
    porofluidmultiphase_ele
    poromultiphase_scatra
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
