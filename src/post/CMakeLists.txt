baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    global_legacy_module
    inpar
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
