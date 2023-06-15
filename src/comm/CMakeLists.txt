baci_auto_define_module()

set(_comm_dependencies
    # cmake-format: sortable
    io
    lib
    linalg
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_comm_dependencies})
