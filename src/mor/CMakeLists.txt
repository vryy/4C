baci_auto_define_module()

set(_mor_dependencies
    # cmake-format: sortable
    adapter
    io
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_mor_dependencies})
