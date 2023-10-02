baci_auto_define_module()

set(LIBS_REBALANCE
    # cmake-format: sortable
    io
    lib
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_REBALANCE})
