baci_auto_define_module()

set(_red_airways_dependencies
    # cmake-format: sortable
    adapter
    core
    inpar
    io
    lib
    mat
    rebalance
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_red_airways_dependencies})
