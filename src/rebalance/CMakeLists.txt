baci_auto_define_module()

set(_rebalance_dependencies
    # cmake-format: sortable
    inpar
    io
    lib
    linalg
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_rebalance_dependencies})
